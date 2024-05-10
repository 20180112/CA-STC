# Optimizing Steganographic Fidelity: Content-Aware Syndrome Trellis Code

**Authors:** Huiyi Tang, Junlong Mao, Shanxiang Lyu, Ling Liu

*College of Cyber Security, Jinan University, Guangzhou 510632, China*

*Guangzhou Institute of Technology, Xidian University, Guangzhou, China*


## PATHON Codes
The PATHON codes compare different STC alorithms including:

- dictionaryCarrier: generate a dictionary for numbering each carrier, with tuple type keys and carrier serial numbers as key values
- messagetravelCarrier: the distortion dictionary obtained by traversing each type of message for each carrier, with the dictionary key being the carrier number and the key value being the distortion dictionary corresponding to each type of message
- status: generate a match dictionary with keys representing the number of each carrier and key values representing the sequence of hidden information exchange corresponding to each carrier
- update_message: update message



## Abstract
Syndrome Trellis Code (STC) stands out as one of the nearly optimal steganographic methods to date. Its superior efficiency and performance have garnered significant attention. However, STC has a limitation: it struggles to handle unevenly distributed host signals and messages, preventing it from achieving the theoretical minimum distortion. To address this flaw, we introduce an enhanced STC scheme called Content-Aware STC (CA-STC). In our work, we propose modifying the codebook based on the statistical relationship between host signals and messages. This adjustment aims to reduce overall distortion. Additionally, we introduce a new parameter to strike a balance between complexity and embedding efficiency in the proposed method. Simulation results, including scenarios involving random data and image datasets, demonstrate that our approach outperforms STC in terms of global distortion, effectively adapting to diverse scenarios.

## Citation
Huiyi Tang, Junlong Mao, Shanxiang Lyu, Ling Liu: "Optimizing Steganographic Fidelity: Content-Aware Syndrome Trellis Code", 2024.
