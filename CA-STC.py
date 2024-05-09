import sys
import random
import numpy as np


class STC:
    def __init__(self, stcode, k):
        """
        初始化
        :param stcode: 十进制形式的子矩阵构造
        :param k: 一次隐藏信息的比特数，也将作为构造校验矩阵大小的参数
        """
        n_bits = len(bin(np.max(stcode))[2:])  # 转为二进制时的长度，即子矩阵的行数h

        submatrix = []
        for d in stcode:
            submatrix.append([int(x) for x in list(bin(d)[2:].ljust(n_bits, '0'))])
        submatrix = np.array(submatrix).T #生成子矩阵
        rows, cols = submatrix.shape

        self.matrix = np.zeros((k + rows - 1, cols * k))  # 用子矩阵构建校验矩阵，可证明其行数和列数一定满足这样的关系
        for i in range(k):
            self.matrix[i: i + rows, cols * i: cols * (i + 1)] = submatrix

        self.matrix_cols = cols * k  # 校验矩阵的列数，也为一次嵌入时载体的长度
        self.submatrix_rows = n_bits  # 子矩阵的行数h
        self.code_h = np.tile(stcode, k)  # 将校验矩阵以列形式记录
        self.code_shift = np.tile([0] * (cols - 1) + [1], k) 
        pass

    def _dual_viterbi(self, x, w, m):
        C = np.zeros((2 ** self.submatrix_rows, self.matrix_cols))
        costs = np.infty * np.ones((2 ** self.submatrix_rows, 1))
        costs[0] = 0
        paths = np.zeros((2 ** self.submatrix_rows, self.matrix_cols))

        m_id = 0  
        y = np.zeros(x.shape)

        for i in range(self.matrix_cols):
            costs_old = costs.copy()
            hi = self.code_h[i]
            ji = 0
            for j in range(2 ** self.submatrix_rows):
                c1 = costs_old[ji] + x[i] * w[i]
                c2 = costs_old[(ji ^ hi)] + (1 - x[i]) * w[i]
                if c1 < c2:
                    costs[j] = c1
                    paths[j, i] = ji  
                else:
                    costs[j] = c2
                    paths[j, i] = ji ^ hi  
                ji = ji + 1

            for j in range(self.code_shift[i]):
                tail = np.infty * np.ones((2 ** (self.submatrix_rows - 1), 1))
                if m[m_id] == 0:
                    costs = np.vstack((costs[::2], tail))
                else:
                    costs = np.vstack((costs[1::2], tail))

                m_id = m_id + 1

            C[:, i] = costs[:, 0]

        ind = np.argmin(costs)
        min_cost = costs[ind, 0]

        m_id -= 1

        for i in range(self.matrix_cols - 1, -1, -1):
            for j in range(self.code_shift[i]):
                ind = 2 * ind + m[m_id, 0]  
                m_id = m_id - 1

            y[i] = paths[ind, i] != ind
            ind = int(paths[ind, i])

        return y.astype('uint8'), min_cost, paths

    def _calc_syndrome(self, x):
        m = np.zeros((np.sum(self.code_shift), 1))
        m_id = 0
        tmp = 0
        for i in range(self.matrix_cols):
            hi = self.code_h[i]
            if x[i] == 1:
                tmp = hi ^ tmp
            for j in range(self.code_shift[i]):
                m[m_id] = tmp % 2
                tmp //= 2
                m_id += 1
        return m.astype('uint8')

    def _bytes_to_bits(self, m):
        bits = []
        for b in m:
            for i in range(8):
                bits.append((b >> i) & 1)
        return bits

    def _bits_to_bytes(self, m):
        enc = bytearray()
        idx = 0
        bitidx = 0
        bitval = 0
        for b in m:
            if bitidx == 8:
                enc.append(bitval)
                bitidx = 0
                bitval = 0
            bitval |= b << bitidx
            bitidx += 1
        if bitidx == 8:
            enc.append(bitval)
        return bytes(enc)
     
    def extract2(self, stego): #返回message为二进制
        y = stego.flatten()
        message = []
        for i in range(0, len(y), self.matrix_cols):
            y_chunk = y[i:i + self.matrix_cols][:, np.newaxis] % 2
            if len(y_chunk) < self.matrix_cols:
                break
            m_chunk = self._calc_syndrome(y_chunk)
            message += m_chunk[:, 0].tolist()
        return message 

    def embed2(self, cover, costs, message):
        #messgage为二进制列表
        shape = cover.shape
        x = cover.flatten()
        w = costs.flatten()
        ml = np.sum(self.code_shift)
        message_bits = np.array(message)

        i = 0
        j = 0
        y = x.copy()
        while True:
            x_chunk = x[i:i + self.matrix_cols][:, np.newaxis] % 2 
            m_chunk = message_bits[j:j + ml][:, np.newaxis]
            w_chunk = w[i:i + self.matrix_cols][:, np.newaxis]
            y_chunk, min_cost, _ = self._dual_viterbi(x_chunk, w_chunk, m_chunk)
            idx = x_chunk[:, 0] != y_chunk[:, 0]
            y[i:i + self.matrix_cols][idx] += 1
            i += self.matrix_cols
            j += ml
            if i + self.matrix_cols > len(x) or j + ml > len(message_bits):
                break
        return np.vstack(y).reshape(shape)
    
    def dictionaryCarrier(self,carrier,k):
        '''
        carrier:生成的随机分布的载体(array)
        k:一次隐藏信息的数量
        return 返回一个对每种载体进行编号的字典，字典键为tuple类型，键值为载体序号
        '''
        dicCarrier = {}
        count = 0
        for i in range(0,len(carrier),2*k):
            temp = tuple(carrier[i:i+2*k].tolist())
            if temp not in dicCarrier.keys():
                dicCarrier[temp] = count
                count = count + 1
        return dicCarrier

    def messagetravelCarrier(self,dicCarrier,k):
        '''
        dicCarrier:载体编号字典
        k:一次隐藏的信息数量
        return 针对每种载体遍历每种信息得到的失真字典，字典键为载体编号，键值为每种信息对应的失真字典
        '''
        dicCarrierMessage = {} #key为每种载体，值为载体对应的信息失真字典
        #对每种载体进行遍历得到每个信息对应的失真
        for i in dicCarrier.keys(): 
            dicMessageDistortion = {}
            for j in range(0,int(pow(2,k))):
                temp = bin(j)[2:]
                while len(temp) < k:
                    temp = "0" + temp
                m = []
                for p in temp:
                    m.append(int(p))
                m = np.asarray(m)
                #对每个Carrier遍历每种隐藏信息并记录下失真作为值
                x_chunk = np.asarray(i)[:, np.newaxis] % 2
                m_chunk = m[:, np.newaxis]
                w_chunk = np.asarray([0]*(2*k))
                w_chunk = w_chunk[:, np.newaxis]
                y_chunk, min_cost, _ = self._dual_viterbi(x_chunk, w_chunk, m_chunk)
                idx = x_chunk[:, 0] != y_chunk[:, 0]   
                dicMessageDistortion[j] = sum(idx)
            dicCarrierMessage[dicCarrier[i]] = dicMessageDistortion         
        return dicCarrierMessage

    def status(self,carrier,dicCarrier,message,dicCarrierMessage,k):
        '''
        carrier:随机生成的载体
        dicCarrier:载体编号字典
        message:随机生成的信息
        dicCarrierMessage:生成的每种载体对应每种信息失真的字典
        k:一次隐藏的信息数量
        return 一个order字典，键为每种载体的编号，键值为每种载体对应的隐藏信息交换序列
        '''
        #得到每个载体对应的order
        c = carrier.flatten()
        message_bits = message.flatten()
        num = np.zeros((len(dicCarrier),int(pow(2,k)))) #记录数量
        weight = np.zeros((len(dicCarrier),int(pow(2,k))))

        i = 0
        j = 0
        while True:
            #遍历载体和隐藏信息，计算出每种载体对应的隐藏信息的数量和失真
            c_chunk_tuple = tuple(c[i:i+2*k].tolist())
            c_chunk = c[i:i+2*k][:,np.newaxis] % 2
            m_chunk = message_bits[j:j+k][:,np.newaxis]
            m_list = message_bits[j:j+k].tolist()
            m_bin = 0
            for p in m_list:
                m_bin = m_bin*2 + p
            w_chunk = np.asarray([0]*(2*k))
            w_chunk = w_chunk[:, np.newaxis]
            num[dicCarrier[c_chunk_tuple]][m_bin-1] = num[dicCarrier[c_chunk_tuple]][m_bin-1] + 1
            y_chunk, min_cost, _ = self._dual_viterbi(c_chunk, w_chunk, m_chunk)
            idx = c_chunk[:, 0] != y_chunk[:, 0]
            weight[dicCarrier[c_chunk_tuple]][m_bin-1] = sum(idx)

            i+=2*k
            j+=k
            if i+2*k > len(c) or j+k > len(message_bits):
                break
        import math
        #得到总失真矩阵
        total_weight = np.multiply(num,weight)
        dicCarrierMessageOrder = {}
        for i in range(len(total_weight)):
            #遍历每个载体，得到每个载体的order
            ordertemp = []
            temp = dicCarrierMessage[i] #此载体中每种m对应的失真的字典
            for j in range(len(total_weight[0])):
                min_key = min(temp,key=lambda k:temp[k])
                max_distortion_m = np.argmax(total_weight[i,:])
                temp[min_key] = math.inf
                total_weight[i,max_distortion_m] = -math.inf
                ordertemp.append((max_distortion_m,min_key))
            dicCarrierMessageOrder[i] = ordertemp
        return dicCarrierMessageOrder

    def status2(self,carrier,message,k):
        '''
        carrier:随机生成的载体
        message:随机生成的信息
        k:一次隐藏的信息数量
        return 未优化前失真，初始得到的y序列，一个order表示y的修改序列
        '''
        #得到每个载体对应的order
        c = carrier.flatten()
        message_bits = message.flatten()
        y = np.asarray([])

        i = 0
        j = 0
        while True:
            #遍历载体和隐藏信息，计算出每种载体对应的隐藏信息的数量和失真
            c_chunk_tuple = tuple(c[i:i+2*k].tolist())
            c_chunk = c[i:i+2*k][:,np.newaxis] % 2
            m_chunk = message_bits[j:j+k][:,np.newaxis]
            m_list = message_bits[j:j+k].tolist()
            m_bin = 0
            for p in m_list:
                m_bin = m_bin*2 + p
            w_chunk = np.asarray([0]*(2*k))
            w_chunk = w_chunk[:, np.newaxis]
            y_chunk, min_cost, _ = self._dual_viterbi(c_chunk, w_chunk, m_chunk)
            y = np.concatenate((y,y_chunk.flatten()))

            i+=2*k
            j+=k
            if i+2*k > len(c) or j+k > len(message_bits):
                break
        idx = (c[:]%2) != y[:]
        
        from scipy.optimize import linear_sum_assignment
        def matrix_max_value(matrix: np.array):
            matrix = -matrix
            row_ind, col_ind = linear_sum_assignment(matrix)
            max_task = [(i, c) for i, c in enumerate(col_ind)]
            return max_task

        matrix = np.zeros((int(pow(2,2)),int(pow(2,2))))
        for i in range(0,len(y),2): #两个比特为一组变
            c_list = (c[i:i+2] % 2).tolist() 
            y_list = (y[i:i+2] % 2).tolist() 
            c_bin = 0
            y_bin = 0
            for p in c_list:
                c_bin = c_bin*2+p
            for p in y_list:
                y_bin = y_bin*2+p
            matrix[int(y_bin-1)][int(c_bin-1)] = matrix[int(y_bin-1)][int(c_bin-1)] + 1
        order = matrix_max_value(matrix)
        print(matrix)
        print(order)
        return sum(idx),y,order

    def update_y(self,carrier,stego,order):
        stego_change = []
        for i in range(0,len(stego),2):
            stego_list = (stego[i:i+2] % 2).tolist() 
            stego_bin = 0
            for p in stego_list:
                stego_bin = stego_bin*2+p
            temp = ""
            for j in order:
                if stego_bin == j[0]:
                    temp = bin(j[1])[2:]
                    break
            while len(temp) < 2:
                temp = "0"+temp
            for k in temp:
                stego_change.append(int(k))
        stego_change = np.asarray(stego_change)
        idx =  (carrier[:]%2) != stego_change[:]
        #print(sum(idx))
        return sum(idx),stego_change

    def update_message(self,dicCarrier,carrier,dicCarrierMessageOrder,message,k):
        '''
        carrier:随机生成的载体
        dicCarrier:载体编号字典
        message:随机生成的信息
        dicCarrierMessageOrder:每种载体对应的隐藏信息交换序列
        k:一次隐藏的信息数量
        return 更新后的message序列
        '''
        #根据每个载体修改message
        update_m = []
        message_bits = message.flatten()
        
        i = 0
        j = 0
        while True:
            order = dicCarrierMessageOrder[dicCarrier[tuple(carrier[i:i+2*k].tolist())]]
            s = message_bits[j:j+k]
            s_str = ""
            for p in s:
                s_str = s_str + str(p)
            s = int(s_str,2)
            for p in range(len(order)):
                if s == order[p][0]:
                    update_m.append(order[p][1])
                    pass
            i+=2*k
            j+=k
            if i+2*k > len(carrier) or j+k > len(message_bits):
                break
        
        update_m_bin = []
        for i in update_m:
            temp = bin(i)[2:]
            while len(temp) < k:
                temp = "0" + temp
            for j in temp:
                update_m_bin.append(int(j))
        return update_m_bin

class Datagenerate:
    def generateMessage(k,percentage,num):
        message = [] 
        for i in range(num*k):
            r = random.random()
            if r < percentage:
                message.append(0)
            else:
                message.append(1)
            pass
        return np.asarray(message)

    def generateCarrier(k,num,mu,sigma):
        carrier = []
        i = 0
        while i < num*2*k:
            cc = round(random.normalvariate(mu,sigma))
            if -256 <= cc <= 256:
                carrier.append(cc)
                i+=1
                pass
            pass
        return np.asarray(carrier)

# Example
if __name__ == "__main__":

    up_distortion = 0
    distortion = 0          
 
    for i in range(3):
        stcode = [109, 71]  #the submatrix 
        k1 = 4 # the dimension of STC
        k = 4 
        num = 30   
        stc = STC(stcode, k1)
        Date = Datagenerate
        #Step 1:Generate random carriers and message
        carrier = Date.generateCarrier(k,num,0,1)
        messagegener = Date.generateMessage(k,0.9,num)
        #Step 2:the dictionary of carrier
        dicCarrier = stc.dictionaryCarrier(carrier,k)           
        #Step 3:the dictionary of message distortion for each carrier
        dicCarrierMessage = stc.messagetravelCarrier(dicCarrier,k)   
        #Step 4:the dictionary of match for each carrier
        dicCarrierMessageOrder = stc.status(carrier,dicCarrier,messagegener,dicCarrierMessage,k)
        #step 5:update message
        update_mess = stc.update_message(dicCarrier,carrier,dicCarrierMessageOrder,messagegener,k)
        
        costs = np.asarray([0]*(2*k*num)) 
        #STC
        stego = stc.embed2(carrier, costs, messagegener)
        extracted_message = stc.extract2(stego)
        print("distortion:",sum(stego - carrier))
        distortion = distortion + sum(stego - carrier)
        #CA-STC
        stego_update = stc.embed2(carrier, costs, update_mess)
        extracted_message_update = stc.extract2(stego_update)
        print("update distortion:",sum(stego_update - carrier))
        up_distortion = up_distortion + sum(stego_update - carrier)
    
    print("distortion:",distortion/3)
    print("update distortion:",up_distortion/3)
