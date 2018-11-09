
from SloppyCell.ReactionNetworks import *

# give a name
expt = Experiment('LitterMCPA')

# read in data from ML treatment
data = {'base': {
                 'cBtot': {0: (0.0253, 0.0023),
                           4.89: (0.0483, 0.0078),
                           7.84: (0.0402, 0.0116),
                           10.03: (0.0436, 0.0071),
                           13.93: (0.0433, 0.0060),
                           22.81: (0.0705, 0.0173),
                           },
                 'x_2': {0: (0.000133, 0.000014),
                         4.89: (0.000242, 0.000085),
                         7.84: (0.000263, 0.000118),
                         10.03: (0.000263, 0.000021),
                         13.93: (0.000372, 0.000061),
                         22.81: (0.000679, 0.000034),
                         },
                'x_3': {0: (0.0160, 0.0024),
                         4.89: (0.0174, 0.0094),
                         7.84: (0.0373, 0.0093),
                         10.03: (0.0332, 0.0137),
                         13.93: (0.0384, 0.0177),
                         22.81: (0.0565, 0.0180),
                        },
                'cmic': {0: (0.041433, 0.016),
                         4.89: (0.0802, 0.0300),
                         7.84: (0.0708, 0.0238),
                         10.03: (0.1220, 0.0576),
                         13.93: (0.1266, 0.0244),
                         22.81: (0.1192, 0.0322),
                        },
                'doc': {0: (0.042, 0.006),
                        4.89: (0.1178, 0.036),
                        7.84: (0.0874, 0.0222),
                        10.03: (0.0958, 0.0228),
                        13.93: (0.1072, 0.0106),
                        22.81: (0.1240, 0.0186),
                        },
                'toc': {0: (14.90, 0.81),
                        4.89: (15.3460, 1.0128),
                        7.84: (15.0820, 0.3118),
                        10.03: (16.1080, 0.3814),
                        13.93: (15.5500, 0.4622),
                        22.81: (15.6760, 0.7199),
                        },
                'x_8': {0: (0.00002, 2.35e-06),
                         4.89: (0.00000675, 1.3e-06),
                         7.84: (0.00000134, 8.2e-7),
                         10.03: (0.0000004477, 3.3000e-07),
                         13.93: (3.7000e-08, 1.0000e-08),
                         22.81: (3.7000e-08, 1.0000e-08)
                         },
                'co2_tot': {
                         0.83: (0.0270, 2.3476e-02),
                         1.79: (0.0330, 0.014),
                         3.76: (0.1144, 0.0191),
                         4.86: (0.1606, 1.9614e-02),
                         5.77: (0.1933, 1.0355e-02),
                         7.74: (0.2599, 1.2809e-02),
                         8.75: (0.2927, 0.0191),
                         9.97: (0.3308, 0.0195),
                         13.76: (0.4239, 0.0132),
                         15.73: (0.4668, 1.8987e-02),
                        18.77: (0.5257, 0.0194),
                        20.76: (0.5606, 0.019),
                         22.74: (0.5718, 1.5361e-02)
                         }
                }
        }

# attach data dictionary
expt.set_data(data)

# no scaling factors
expt.set_fixed_sf({'x_1':1, 'x_2':1, 'x_3':1,'x_4':1, 'x_5':1, 'x_6':1, 'x_7':1, 'x_8':1, 'x_9':1, 'x_10':1,'x_11':1, 'x_12':1,
                   'cmic':1, 'cBtot':1, 'doc':1, 'toc':1, 'co2_tot':1})
