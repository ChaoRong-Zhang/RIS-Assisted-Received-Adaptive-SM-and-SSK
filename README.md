# RIS-Assisted-RASM-RASSK

This repository contains MATLAB simulation code for RIS-assisted received adaptive spatial modulation (RASM) and received adaptive space shift keying (RASSK), associated with our IEEE TVT 2026 and IEEE WCNC 2025 papers.

## 📖 Overview

The proposed RASM and RASSK schemes incorporate adaptive receive-antenna selection into a reconfigurable intelligent surface (RIS)-assisted wireless communication system. Information bits are mapped to receive-antenna combinations, with the number of selected antennas varying across transmitted symbols.

RASM conveys information through both antenna-combination indices and constellation symbols, whereas RASSK conveys information through antenna-combination indices without additional constellation modulation. The RIS adjusts its reflection phases according to the selected antenna combination.

The TVT paper extends the earlier WCNC work by introducing RASSK and further analyzing reliability and physical-layer security. This code release includes spectral-efficiency calculations and Monte Carlo bit-error-rate (BER) simulations for the proposed schemes and comparison baselines.

## 📄 Papers

**IEEE TVT 2026 — Journal Paper**

From Reliability to Security: How RIS-Assisted Adaptive SM and SSK Enhances Wireless Systems

[IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/11506249) | [arXiv](https://arxiv.org/abs/2512.03518)

**IEEE WCNC 2025 — Conference Paper**

RIS-Assisted Received Adaptive Spatial Modulation for Wireless Communications

[IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/10978803) | [arXiv](https://arxiv.org/abs/2407.06894)

## 📝 Citation

If you find this repository helpful for your academic research, please consider citing the relevant paper(s):

```bibtex
@article{zhang2026reliability,
  title={From Reliability to Security: How {RIS}-Assisted Adaptive {SM} and {SSK} Enhances Wireless Systems},
  author={Zhang, Chaorong and Ng, Benjamin K. and Wang, Ke and Xu, Hui and Lam, Chan-Tong},
  journal={IEEE Transactions on Vehicular Technology},
  year={2026},
  publisher={IEEE},
  doi={10.1109/TVT.2026.3690516}
}

@inproceedings{zhang2025ris,
  title={{RIS}-Assisted Received Adaptive Spatial Modulation for Wireless Communications},
  author={Zhang, Chaorong and Xu, Hui and Ng, Benjamin K. and Lam, Chan-Tong and Wang, Ke},
  booktitle={2025 IEEE Wireless Communications and Networking Conference (WCNC)},
  year={2025},
  publisher={IEEE},
  doi={10.1109/WCNC61545.2025.10978803}
}
```

**IEEE format:**

C. Zhang, B. K. Ng, K. Wang, H. Xu, and C.-T. Lam, “From reliability to security: How RIS-assisted adaptive SM and SSK enhances wireless systems,” *IEEE Transactions on Vehicular Technology*, 2026, doi: 10.1109/TVT.2026.3690516.

C. Zhang, H. Xu, B. K. Ng, C.-T. Lam, and K. Wang, “RIS-assisted received adaptive spatial modulation for wireless communications,” in *2025 IEEE Wireless Communications and Networking Conference (WCNC)*, 2025, doi: 10.1109/WCNC61545.2025.10978803.

## 📜 License

This code is provided for non-commercial academic research under the terms in `LICENSE`. Non-commercial academic use, modification, and redistribution are permitted subject to those terms.

Commercial use requires separate written permission from the copyright holders. The code license does not apply to the associated papers or publisher-owned materials.

## 📜 Email Address
zcryyds666@gmail.com
