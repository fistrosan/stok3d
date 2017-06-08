MODULE DAMPING
  !
  ! J M Borrero
  ! July 22, 2013
  ! KIS, Freiburg
  !
  USE CONS_PARAM
  USE ATOM_DATABASE
  USE LOG
  USE SPLINES
  USE DERIVVAR
  !  
  !  This is the s-p data from Anstee and O'Mara 1995, MNRAS 276,859
  !  
  REAL(SP), PARAMETER :: CS_SP(18,21) = RESHAPE(SHAPE= (/18,21 /), SOURCE= &
       (/ 126,   140,   165,  202,  247,  299,  346,  383,  435,  491,  553,  617,  685,  769,  838,  925, 1011, 1082,&
       140,   150,   162,  183,  218,  273,  327,  385,  440,  501,  557,  620,  701,  764,  838,  923, 1025, 1085,&
       154,   167,   175,  192,  216,  251,  299,  357,  423,  487,  549,  617,  684,  759,  834,  910, 1014, 1064,&
       166,   180,   192,  206,  226,  253,  291,  339,  397,  459,  532,  600,  676,  755,  832,  896, 1002, 1055,&
       208,   194,   207,  223,  242,  265,  296,  335,  384,  445,  511,  583,  656,  726,  817,  889,  988, 1044,&
       262,   254,   220,  239,  261,  283,  310,  344,  388,  442,  496,  568,  635,  725,  791,  890,  970, 1036,&
       311,   306,   299,  251,  280,  304,  330,  361,  396,  443,  500,  563,  630,  704,  796,  880,  951, 1033,&
       358,   359,   350,  338,  293,  323,  352,  381,  416,  455,  511,  566,  635,  706,  780,  859,  946, 1039,&
       411,   409,   405,  392,  370,  340,  375,  406,  439,  478,  525,  580,  644,  714,  790,  873,  961, 1050,&
       462,   463,   459,  450,  443,  400,  394,  432,  467,  501,  546,  595,  650,  711,  786,  873,  963, 1050,&
       522,   525,   529,  524,  516,  518,  438,  454,  495,  532,  565,  621,  671,  741,  813,  874,  951, 1034,&
       589,   593,   590,  583,  579,  568,  565,  483,  517,  560,  600,  644,  691,  752,  821,  904,  978, 1048,&
       658,   655,   666,  657,  649,  653,  649,  587,  549,  592,  674,  674,  728,  782,  833,  902,  992, 1084,&
       738,   742,   747,  725,  721,  729,  699,  730,  626,  622,  668,  721,  765,  809,  887,  938, 1001, 1109,&
       838,   838,   810,  809,  790,  800,  769,  815,  757,  679,  704,  755,  806,  854,  901,  974, 1034, 1105,&
       942,   946,   925,  901,  918,  895,  919,  897,  933,  890,  785,  797,  859,  908,  976, 1020, 1115, 1173,&
       1059,  1061,  1056, 1061, 1074, 1031, 1036, 1036,  993, 1038,  932,  852,  878,  943, 1003, 1074, 1131, 1200,&
       1069,  1076,  1083, 1095, 1102, 1091, 1126, 1156, 1103, 1149, 1157, 1036,  972, 1007, 1064, 1124, 1209, 1283,&
       1338,  1350,  1356, 1354, 1324, 1301, 1312, 1318, 1257, 1239, 1297, 1233, 1089, 1059, 1106, 1180, 1218, 1317,&
       1409,  1398,  1367, 1336, 1313, 1313, 1409, 1354, 1317, 1287, 1353, 1386, 1279, 1158, 1141, 1188, 1260, 1335,&
       1328,  1332,  1342, 1369, 1405, 1451, 1502, 1524, 1506, 1477, 1522, 1594, 1572, 1436, 1328, 1325, 1382, 1446 /))
  !
  REAL(SP), PARAMETER :: AL_SP(18,21) = RESHAPE(SHAPE= (/18,21 /), SOURCE= &
       (/.268, .269, .335, .377, .327, .286, .273, .270, .271, .268, .267, .264, .264, .264, .261, .256, .248, .245,&
       .261, .256, .254, .282, .327, .355, .321, .293, .287, .271, .267, .273, .270, .270, .268, .268, .264, .263,&
       .266, .264, .257, .252, .267, .289, .325, .339, .319, .301, .292, .284, .281, .281, .277, .282, .276, .274,&
       .262, .274, .258, .251, .247, .254, .273, .291, .316, .322, .320, .302, .294, .290, .287, .292, .283, .277,& 
       .322, .275, .264, .259, .250, .245, .273, .255, .271, .284, .294, .308, .296, .299, .288, .289, .282, .278,&
       .267, .300, .260, .268, .254, .242, .243, .242, .239, .246, .267, .277, .280, .290, .282, .281, .274, .271,&
       .259, .274, .275, .252, .265, .248, .249, .237, .283, .236, .247, .254, .254, .271, .268, .267, .258, .262,&
       .260, .255, .268, .268, .268, .264, .248, .239, .229, .240, .236, .234, .238, .244, .252, .251, .244, .255,&
       .255, .255, .244, .247, .317, .246, .255, .244, .237, .231, .227, .231, .235, .232, .235, .241, .237, .245,&
       .256, .254, .254, .249, .227, .319, .253, .253, .240, .237, .238, .233, .231, .230, .228, .234, .227, .241,&
       .257, .254, .252, .235, .253, .240, .284, .251, .246, .241, .235, .228, .222, .225, .225, .219, .228, .233,&
       .244, .240, .245, .238, .248, .230, .283, .252, .244, .244, .238, .235, .234, .236, .228, .224, .225, .231,&
       .244, .241, .244, .237, .237, .249, .219, .324, .239, .245, .242, .242, .232, .233, .221, .227, .231, .218,&
       .241, .245, .249, .239, .243, .250, .217, .254, .308, .237, .247, .244, .234, .228, .233, .224, .227, .226,&
       .243, .243, .232, .227, .235, .253, .227, .220, .320, .270, .243, .252, .248, .238, .234, .241, .225, .227,&
       .225, .226, .234, .230, .226, .233, .249, .225, .216, .300, .286, .237, .240, .247, .243, .234, .231, .238,&
       .268, .260, .247, .238, .233, .241, .254, .248, .207, .227, .315, .260, .226, .237, .240, .239, .239, .240,&
       .248, .246, .238, .226, .213, .221, .226, .226, .204, .194, .248, .316, .234, .216, .236, .233, .221, .230,&
       .200, .202, .198, .194, .206, .207, .227, .224, .207, .185, .198, .275, .315, .233, .229, .231, .233, .236,&
       .202, .209, .221, .226, .230, .245, .202, .257, .246, .225, .215, .246, .320, .321, .244, .239, .251, .253,& 
       .246, .248, .255, .265, .274, .285, .292, .284, .273, .250, .225, .239, .295, .352, .320, .258, .260, .269/))   
  ! 
  !  p-d data from Barklem and O'Mara 1997, MNRAS, 290, 102
  !
  REAL(SP), PARAMETER :: CS_PD(18,18) =  RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/425,  461,  507,  566,  630,  706,  799,  889,  995, 1083, 1191, 1334, 1478, 1608, 1790, 1870, 1936, 2140,&
       429,  460,  505,  565,  633,  704,  795,  896,  985, 1082, 1199, 1340, 1487, 1611, 1795, 1872, 1937, 2136,&
       419,  451,  501,  556,  627,  700,  785,  891,  977, 1088, 1212, 1346, 1493, 1604, 1793, 1863, 1930, 2144,&
       402,  437,  489,  544,  614,  695,  779,  875,  975, 1102, 1221, 1350, 1488, 1591, 1774, 1844, 1919, 2126,&
       384,  418,  467,  529,  595,  674,  769,  856,  976, 1108, 1224, 1338, 1467, 1570, 1743, 1817, 1900, 2118,&
       366,  397,  443,  505,  576,  651,  755,  841,  973, 1095, 1210, 1308, 1435, 1545, 1702, 1786, 1878, 2081,&
       356,  387,  432,  489,  562,  635,  722,  841,  961, 1078, 1175, 1273, 1397, 1517, 1672, 1763, 1863, 2034,&
       359,  388,  431,  479,  545,  624,  707,  834,  943, 1059, 1158, 1256, 1368, 1490, 1647, 1747, 1849, 1998,&
       361,  394,  436,  483,  547,  615,  704,  817,  920, 1027, 1124, 1238, 1358, 1465, 1624, 1736, 1838, 1978,&
       400,  382,  440,  489,  546,  610,  690,  817,  897,  998, 1115, 1201, 1351, 1453, 1599, 1728, 1829, 1953,&
       474,  461,  416,  491,  549,  612,  701,  806,  883,  974, 1078, 1194, 1310, 1456, 1569, 1716, 1818, 1925,&
       531,  518,  507,  463,  547,  615,  694,  784,  881,  958, 1047, 1153, 1297, 1432, 1547, 1688, 1809, 1901,&
       594,  585,  577,  564,  513,  615,  695,  779,  879,  949, 1041, 1145, 1264, 1388, 1544, 1644, 1804, 1879,&
       675,  659,  651,  639,  632,  576,  695,  782,  879,  957, 1046, 1141, 1254, 1391, 1524, 1614, 1793, 1871,&
       739,  734,  726,  719,  715,  708,  663,  776,  901,  971, 1022, 1117, 1232, 1355, 1478, 1616, 1766, 1887,&
       819,  821,  805,  784,  773,  761,  736,  761,  888,  958, 1044, 1145, 1237, 1346, 1487, 1614, 1721, 1891,&
       899,  895,  871,  852,  856,  861,  854,  759,  883,  984, 1027, 1113, 1226, 1355, 1467, 1568, 1703, 1885,&
       973,  946,  955,  925,  939,  927,  902,  920,  870,  987, 1061, 1145, 1234, 1319, 1439, 1552, 1722, 1859/))
  !
  REAL(SP), PARAMETER :: AL_PD(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/.281, .288, .283, .282, .278, .281, .272, .274, .268, .257, .251, .243, .246, .251, .254, .268, .304, .308,&  
       .290, .297, .291, .290, .286, .282, .277, .275, .267, .254, .252, .244, .250, .257, .260, .274, .308, .312,& 
       .294, .299, .293, .294, .288, .289, .281, .276, .265, .256, .251, .247, .258, .264, .268, .283, .318, .317,&  
       .297, .298, .302, .300, .289, .295, .290, .276, .264, .256, .260, .258, .268, .277, .281, .292, .330, .327,&  
       .305, .311, .313, .315, .305, .304, .299, .279, .271, .272, .273, .276, .285, .290, .293, .302, .340, .340,&  
       .292, .294, .303, .305, .301, .307, .290, .277, .274, .278, .287, .288, .295, .302, .306, .312, .343, .346,&  
       .268, .277, .279, .285, .285, .290, .279, .278, .280, .283, .295, .296, .305, .310, .313, .315, .342, .346,&  
       .288, .285, .280, .278, .278, .277, .272, .271, .279, .288, .297, .305, .310, .313, .311, .310, .335, .338,&  
       .314, .304, .292, .282, .275, .275, .262, .272, .290, .293, .299, .307, .308, .310, .303, .302, .325, .328,&  
       .346, .329, .313, .295, .283, .275, .264, .274, .288, .302, .307, .310, .306, .307, .292, .296, .315, .320,&  
       .320, .295, .326, .318, .294, .277, .275, .271, .293, .303, .305, .309, .309, .303, .294, .294, .310, .313,&  
       .304, .310, .297, .320, .317, .297, .283, .274, .298, .305, .308, .311, .313, .300, .290, .293, .305, .306,&  
       .314, .313, .308, .297, .325, .314, .293, .276, .292, .309, .314, .308, .303, .296, .286, .291, .301, .302,&  
       .308, .311, .307, .312, .288, .340, .305, .285, .294, .310, .315, .309, .296, .285, .281, .288, .298, .295,&  
       .313, .310, .315, .303, .313, .294, .331, .286, .294, .307, .320, .316, .303, .281, .278, .285, .290, .292,&  
       .315, .306, .308, .297, .295, .283, .334, .297, .280, .294, .314, .321, .313, .291, .280, .279, .287, .290,&  
       .308, .304, .305, .297, .279, .285, .251, .278, .278, .284, .297, .314, .307, .289, .274, .274, .274, .291,& 
       .301, .299, .298, .285, .265, .279, .241, .285, .260, .286, .302, .306, .302, .288, .277, .263, .271, .293/))
  
  !  
  ! d-f data from Barklem, O'Mara and Ross, 1998, MNRAS, 296, 1057
  !
  REAL(SP), PARAMETER :: CS_DF(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/808,  873,  958, 1059, 1175, 1306, 1453, 1615, 1793, 1979, 2121, 2203, 2461, 2604, 2764, 2757, 2784, 3156,&
       798,  866,  953, 1052, 1172, 1299, 1450, 1606, 1776, 1967, 2114, 2196, 2451, 2601, 2763, 2767, 2783, 3142,&
       781,  848,  934, 1030, 1149, 1276, 1416, 1596, 1751, 1944, 2100, 2188, 2436, 2594, 2767, 2777, 2795, 3123,&
       766,  831,  915, 1010, 1124, 1239, 1398, 1564, 1729, 1912, 2083, 2180, 2426, 2585, 2776, 2790, 2808, 3106,&
       750,  814,  897,  987, 1097, 1201, 1355, 1530, 1718, 1875, 2060, 2171, 2414, 2575, 2779, 2809, 2820, 3103,&
       733,  797,  872,  950, 1049, 1166, 1326, 1502, 1670, 1851, 2026, 2165, 2396, 2562, 2779, 2827, 2832, 3099,&
       726,  786,  853,  936, 1011, 1128, 1303, 1472, 1649, 1844, 1979, 2159, 2371, 2548, 2778, 2840, 2848, 3103,&
       709,  783,  847,  912, 1002, 1093, 1270, 1419, 1606, 1787, 1951, 2139, 2335, 2533, 2775, 2847, 2863, 3104,&
       758,  721,  838,  907, 1010, 1066, 1211, 1401, 1600, 1774, 1972, 2098, 2313, 2528, 2781, 2857, 2892, 3121,&
       869,  882,  820,  870, 1003, 1098, 1165, 1368, 1527, 1735, 1896, 2030, 2288, 2534, 2776, 2844, 2902, 3123,&
       970,  967,  934,  938,  918, 1130, 1194, 1287, 1507, 1679, 1821, 2021, 2271, 2525, 2732, 2786, 2882, 3085,&
       1079, 1043, 1056, 1007, 1014, 1021, 1200, 1326, 1424, 1668, 1818, 1988, 2242, 2493, 2672, 2719, 2853, 3035,&
       1174, 1173, 1127, 1154, 1104, 1099, 1169, 1288, 1442, 1580, 1704, 1882, 2136, 2400, 2561, 2648, 2832, 2994,&
       1285, 1278, 1269, 1225, 1252, 1229, 1116, 1343, 1380, 1594, 1710, 1874, 2054, 2309, 2484, 2607, 2813, 2932,&
       1440, 1408, 1422, 1380, 1383, 1341, 1361, 1192, 1448, 1454, 1675, 1873, 2069, 2246, 2432, 2610, 2811, 2878,&
       1572, 1545, 1553, 1517, 1481, 1502, 1469, 1349, 1373, 1561, 1586, 1781, 2072, 2301, 2490, 2626, 2754, 2832,&
       1698, 1701, 1694, 1641, 1617, 1651, 1566, 1600, 1374, 1547, 1698, 1749, 1989, 2289, 2511, 2594, 2689, 2774,&
       1870, 1841, 1786, 1752, 1777, 1757, 1666, 1732, 1522, 1533, 1707, 1817, 1928, 2194, 2435, 2574, 2665, 2742/))
  !
  REAL(SP), PARAMETER :: AL_DF(18,18) = RESHAPE(SHAPE= (/18,18 /), SOURCE= &
       (/.295, .286, .299, .300, .307, .310, .311, .311, .316, .319, .325, .351, .364, .369, .372, .379, .373, .351,&
       .295, .295, .301, .302, .311, .316, .314, .314, .320, .321, .324, .349, .361, .365, .368, .374, .368, .349,&  
       .286, .298, .302, .304, .311, .323, .321, .319, .324, .323, .323, .345, .355, .358, .362, .367, .361, .343,& 
       .290, .295, .307, .316, .322, .329, .326, .325, .329, .324, .321, .343, .350, .351, .354, .360, .358, .337,&  
       .292, .299, .307, .321, .327, .336, .333, .330, .330, .320, .321, .338, .344, .344, .345, .352, .352, .332,&  
       .291, .299, .309, .323, .335, .339, .335, .333, .327, .323, .319, .333, .336, .336, .336, .344, .345, .329,& 
       .297, .302, .312, .321, .340, .338, .333, .327, .325, .319, .318, .324, .329, .330, .330, .336, .337, .325,& 
       .319, .314, .317, .327, .334, .344, .339, .327, .323, .318, .312, .318, .319, .322, .322, .326, .327, .316,&  
       .333, .328, .339, .325, .359, .351, .332, .325, .322, .311, .309, .310, .311, .316, .314, .317, .321, .313,& 
       .274, .273, .323, .412, .318, .339, .359, .328, .324, .311, .309, .325, .322, .315, .318, .319, .325, .314,&  
       .297, .296, .273, .302, .436, .325, .354, .335, .326, .311, .314, .330, .323, .324, .325, .323, .330, .314,& 
       .284, .295, .296, .280, .300, .438, .322, .348, .332, .318, .320, .332, .335, .334, .335, .331, .333, .309,& 
       .280, .278, .285, .297, .279, .320, .445, .319, .320, .324, .328, .338, .348, .346, .345, .336, .328, .300,& 
       .280, .273, .267, .273, .284, .268, .343, .390, .323, .308, .318, .325, .343, .348, .346, .337, .311, .286,&  
       .277, .270, .260, .266, .276, .263, .294, .408, .337, .324, .299, .308, .331, .334, .345, .327, .315, .280,& 
       .270, .262, .258, .260, .273, .273, .262, .375, .410, .298, .312, .294, .313, .331, .328, .322, .307, .270,&  
       .271, .267, .262, .264, .274, .269, .261, .323, .351, .359, .294, .325, .310, .318, .321, .315, .291, .268,& 
       .275, .276, .272, .276, .279, .270, .264, .295, .393, .340, .319, .287, .320, .330, .316, .302, .280, .261 /))
  !
CONTAINS
  ! get_damping
  ! get_neff
  ! get_collisional_param
  ! transition_sp
  ! transition_pd
  ! transition_df
  ! radiative_damping
  !
  ! ---------------------------
  ! Subroutine GET_DAMPING
  ! ----------------------------
  SUBROUTINE GET_DAMPING(I,LAMBDA,ION,EPLOW,NHYD,TEMP,SIGMA,ALPHA,DAMP,DADT_P,DADT_R)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)              :: I, ION
    REAL(SP),   INTENT(IN)              :: TEMP
    REAL(DP),   INTENT(IN)              :: EPLOW, NHYD, SIGMA, ALPHA, LAMBDA
    REAL(SP),   INTENT(OUT)             :: DAMP
    REAL(DP),   INTENT(OUT), OPTIONAL   :: DADT_P, DADT_R
    ! Internal
    REAL(SP),   PARAMETER               :: RAT_POL_HE_H = 0.3181818
    REAL(SP),   PARAMETER               :: RAT_POL_H2_H = 1.2121212
    REAL(SP),   PARAMETER               :: V0=1E6 ! cm/s
    REAL(SP)                            :: MU_X_H, MU_X_HE, E1, E2, GAMMAF, BETA
    REAL(SP)                            :: GAMMA_COL, GAMMA_RAD, VDOPPLER
    REAL(DP)                            :: NHE
    REAL(DP)                            :: DVDOPDTEMP, DGAMMACDTEMP_PG, DGAMMACDTEMP_RHO
    !
    ! Determine Doppler width (in velocity units: cm/s)
    !
    VDOPPLER = SQRT(2.*REAL(KBOL)*TEMP/(REAL(MAMU)*MATOM(I)*SPI))
    DVDOPDTEMP = DBLE(VDOPPLER)/(2D0*DBLE(TEMP))
    !
    ! Determine radiative damping
    CALL RADIATIVE_DAMPING(LAMBDA,GAMMA_RAD)
    !
    !
    ! Calculate Collisional damping according to ABO Theory: Anstee, Barklem & O'Mara
    !
    IF (SIGMA.GT.0.AND.ALPHA.GT.0) THEN
       ! Reduced mass of atom under consideration and Hydrogen
       MU_X_H = REAL(MAMU)* (MATOM(I)*MATOM(1))/(MATOM(I)+MATOM(1))
       ! Reduced mass of atom under consideration and Helium
       MU_X_HE = REAL(MAMU)* (MATOM(I)*MATOM(2))/(MATOM(I)+MATOM(2))
       !
       E1 = 2.-REAL(ALPHA)/2.
       E2 = (1.-REAL(ALPHA))/2.
       ! Gamma function: 5th order polynomial interpolation over the [1,2] region
       GAMMAF=3.6924374-6.5872093*E1+6.3611214*E1**2.-3.2531344*E1**3.+0.88197419*E1**4.-0.095283853*E1**5.
       ! Beta factor (if left like this GAMMA_COL will be in HWHM)
       BETA=(4./SPI)**(REAL(ALPHA)/2.)*GAMMAF*(V0**REAL(ALPHA))*REAL(SIGMA*RBOHR**2.)*(8.*REAL(KBOL)/SPI)**E2
       ! Now it will be in FWHM
       BETA=2.*BETA
       ! Calculate number helium atoms per cm^3
       NHE = NHYD*10D0**(ABUND(2)-12D0)
       ! Calculate collisional damping: this is already FWHM
       GAMMA_COL = REAL(BETA*TEMP**E2*(NHYD*MU_X_H**(-E2)+RAT_POL_HE_H*NHE*MU_X_HE**(-E2)))
       ! Units of GAMMA_COL are s^(-1)
       ! Now we deal with the derivatives of the collisional damping with respect to T and constant PG and RHO
       DGAMMACDTEMP_PG = DBLE(GAMMA_COL*E2/TEMP)+DBLE(BETA*TEMP**E2*DNHYDDTEMP_PG*(MU_X_H**(-E2)+RAT_POL_HE_H* &
            MU_X_HE**(-E2)*10D0**(ABUND(2)-12D0)))
       DGAMMACDTEMP_RHO = DBLE(GAMMA_COL*E2/TEMP)+DBLE(BETA*TEMP**E2*DNHYDDTEMP_RHO*(MU_X_H**(-E2)+RAT_POL_HE_H* &
            MU_X_HE**(-E2)*10D0**(ABUND(2)-12D0)))
       !
    END IF
    ! Finally add both to obtain total damping
    DAMP = (GAMMA_RAD+GAMMA_COL)*REAL(LAMBDA)/(4.*SPI*VDOPPLER)
    ! Now derivatives of teh total damping with respect to T and constant PG and RHO
    DDAMDTEMP_PG = DBLE(DAMP/(GAMMA_RAD+GAMMA_COL))*DGAMMACDTEMP_PG-DBLE(DAMP/VDOPPLER)*DVDOPDTEMP
    DDAMDTEMP_RHO = DBLE(DAMP/(GAMMA_RAD+GAMMA_COL))*DGAMMACDTEMP_RHO-DBLE(DAMP/VDOPPLER)*DVDOPDTEMP
    !
    IF (PRESENT(DADT_P).AND.PRESENT(DADT_R)) THEN
       DADT_P = DDAMDTEMP_PG
       DADT_R = DDAMDTEMP_RHO
    ENDIF
    !
    !
    ! Calculate Collisional damping according to Unsold theory: following Witmann 1972
    !
    !IF (SIGMA.EQ.0.OR.ALPHA.EQ.0) THEN
    !   IF (ION.EQ.1) THEN 
    !      E1=MAXVAL(XI(I)-EPLOW-XI(I)*REAL(ION-1.),1.)
    !      E2=MAXVAL(XI(I)-EPLOW-XI(I)*REAL(ION-1.),3.)
    !   ENDIF
    !   IF (ION.EQ.2) THEN
    !      E1=MAXVAL(XII(I)-EPLOW-XI(I)*REAL(ION-1.),1.)
    !      E2=MAXVAL(XII(I)-EPLOW-XI(I)*REAL(ION-1.),3.)
    !   ENDIF
    !   CHYD = LAMBDA*10.**(0.4*ALOG10(1./E1**2.-1./E2**2.)-12.213)*5.34784E3
    !   IF (ION.EQ.2) CHYD=CHYD*1.741
    !ENDIF
  END SUBROUTINE GET_DAMPING
  ! ----------------------------
  ! Subroutine GET_NEFF
  ! ----------------------------
  SUBROUTINE GET_NEFF(ZN,EPLOW,L0,ION,NLOW_EFF,NUPP_EFF)
    !
    IMPLICIT NONE
    INTEGER,  INTENT(IN)     :: ZN, ION
    REAL(DP), INTENT(IN)     :: EPLOW, L0
    REAL(DP), INTENT(OUT)    :: NLOW_EFF, NUPP_EFF
    REAL(DP)                 :: XH, XLIMIT, XLOW, XUPP
    !
    XLOW=EPLOW*EVOLT                  ! Energy (erg) of the lower level
    XUPP=XLOW+HPLA*LIGHT/(L0*1E-8)    ! Energy (erg) of the upper level
    XH=XII(1)*EVOLT                   ! Ionization potential of hydrogen
    ! -------------------------------------------------------------------------
    ! Series Limit: not exactly ionization potential, but in general very close
    ! -------------------------------------------------------------------------
    ! Single ionized elements Z > 1
    IF (ION.EQ.1.AND.ZN.GT.1) XLIMIT=DBLE(XI(ZN))*EVOLT
    ! Double ionized elements Z > 1
    IF (ION.EQ.2.AND.ZN.GT.1) XLIMIT=DBLE(XII(ZN))*EVOLT
    ! Case of hydrogen Z = 1
    IF (ZN.EQ.1) XLIMIT=DBLE(XII(ZN))*EVOLT
    IF (ZN.EQ.1.AND.ION.GT.1) THEN
       PRINT*,'Hydrogen cannot be ionized twice !'
       PRINT*,'Error in file lines_database.dat. STOP'
       STOP
    ENDIF
    ! Converting from ergs to cm-1
    XLOW=XLOW/(HPLA*LIGHT)
    XUPP=XUPP/(HPLA*LIGHT)
    XH=XH/(HPLA*LIGHT)
    XLIMIT=XLIMIT/(HPLA*LIGHT)
    !
    ! Effective quantum numbers
    !
    NLOW_EFF=SQRT(XH/(XLIMIT-XLOW))
    NUPP_EFF=SQRT(XH/(XLIMIT-XUPP))
    !
  END SUBROUTINE GET_NEFF
  ! --------------------------------
  ! Subroutine GET_COLLISIONAL_PARAM
  ! --------------------------------
  SUBROUTINE GET_COLLISIONAL_PARAM(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    !
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(INOUT) :: IERR
    INTEGER                     :: TRANSITION_TYPE
    !
    ! IERR values = 1 (not a sp, pd, df transition)
    !             = 2 (lower level outside table ranges)
    !             = 3 (upper level outside table ranges)
    TRANSITION_TYPE = SUM(OTRANSITION)
    SELECT CASE(TRANSITION_TYPE)
    CASE(1)
       ! s-p
       CALL TRANSITION_SP(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE(3)
       ! p-d
       CALL TRANSITION_PD(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE(5)
       ! d-f
       CALL TRANSITION_DF(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    CASE DEFAULT
       IERR = 1
    END SELECT
    !print*,nlow_eff,nupp_eff,otransition,sigma,alpha,ierr
    !
  END SUBROUTINE GET_COLLISIONAL_PARAM
  !-------------------------
  ! Subroutine transition_sp
  !-------------------------
  SUBROUTINE TRANSITION_SP(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    !
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(INOUT) :: IERR
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    !
    REAL(DP)     :: CS(21,18), AL(21,18), CS2(21,18), AL2(21,18)
    REAL(DP)     :: NSSG(21), NSPG(18), NSS, NSP
    INTEGER      :: I, J
    !
    IF (OTRANSITION(1) .EQ. 0) NSS=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 1) NSP=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 0) NSS=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 1) NSP=NUPP_EFF
    ! Check table limits
    IF ((NSS .GT. 3.) .OR. (NSS .LT. 1.) .OR. (NSP .GT. 3.) .OR. (NSP .LT. 1.3)) THEN
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       !
       DO I=1,21
          NSSG(I)=1.0+REAL(I-1)*0.1
       ENDDO
       !
       DO I=1,18
          NSPG(I)=1.3+REAL(I-1)*0.1
       ENDDO
       !
       DO I=1,21
          DO J=1,18
             CS(I,J)=CS_SP(J,I)
             AL(I,J)=AL_SP(J,I)
          ENDDO
       ENDDO
       ! second derivative table for spline
       CALL SPLIE2(NSSG,NSPG,CS,21,18,CS2)
       CALL SPLIE2(NSSG,NSPG,AL,21,18,AL2)
       ! run bicubic spline interpolation
       CALL SPLIN2(NSSG,NSPG,CS,CS2,21,18,NSS,NSP,SIGMA)
       CALL SPLIN2(NSSG,NSPG,AL,AL2,21,18,NSS,NSP,ALPHA)
       !
    ENDIF
    !
  END SUBROUTINE TRANSITION_SP
  ! ------------------------
  ! Subroutine transition_pd
  ! ------------------------
  SUBROUTINE TRANSITION_PD(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA ,ALPHA
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    INTEGER,      INTENT(INOUT) :: IERR
    !
    REAL(DP)     :: CS(18,18), AL(18,18), CS2(18,18), AL2(18,18)
    REAL(DP)     :: NSDG(18), NSPG(18), NSD, NSP
    INTEGER      :: I, J
    !
    IF (OTRANSITION(1) .EQ. 1) NSP=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 2) NSD=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 1) NSP=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 2) NSD=NUPP_EFF
    ! Check table limits
    IF ((NSP .GT. 3.) .OR. (NSP .LT. 1.3) .OR. (NSD .GT. 4.) .OR. (NSD .LT. 2.3)) THEN
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       !
       DO I=1,18
          NSPG(I)=1.3+REAL(I-1)*0.1
          NSDG(I)=2.3+REAL(I-1)*0.1
          DO J=1,18
             CS(I,J)=CS_PD(J,I)
             AL(I,J)=AL_PD(J,I)
          ENDDO
       ENDDO
       ! setup second derivative table for spline
       CALL SPLIE2(NSPG,NSDG,CS,18,18,CS2)
       CALL SPLIE2(NSPG,NSDG,AL,18,18,AL2)
       ! run bicubic spline interpolation
       CALL SPLIN2(NSPG,NSDG,CS,CS2,18,18,NSP,NSD,SIGMA)
       CALL SPLIN2(NSPG,NSDG,AL,AL2,18,18,NSP,NSD,ALPHA)
       !
    ENDIF
    !
  END SUBROUTINE TRANSITION_PD
  !
  ! Subroutine transition_df
  !
  SUBROUTINE TRANSITION_DF(NLOW_EFF,NUPP_EFF,OTRANSITION,SIGMA,ALPHA,IERR)
    IMPLICIT NONE
    REAL(DP),     INTENT(IN)    :: NLOW_EFF, NUPP_EFF
    REAL(DP),     INTENT(OUT)   :: SIGMA, ALPHA
    INTEGER,      INTENT(IN)    :: OTRANSITION(2)
    INTEGER,      INTENT(INOUT) :: IERR
    !
    REAL(DP)     :: CS(18,18), AL(18,18), CS2(18,18), AL2(18,18)
    REAL(DP)     :: NSDG(18), NSFG(18), NSD, NSF
    INTEGER      :: I,J
    !
    IF (OTRANSITION(1) .EQ. 2) NSD=NLOW_EFF
    IF (OTRANSITION(1) .EQ. 3) NSF=NLOW_EFF
    IF (OTRANSITION(2) .EQ. 2) NSD=NUPP_EFF
    IF (OTRANSITION(2) .EQ. 3) NSF=NUPP_EFF
    ! Check table limits
    IF ((NSD .GT. 4.) .OR. (NSD .LT. 2.3) .OR. (NSF .GT. 5.) .OR. (NSF .LT. 3.3)) THEN
       IERR=3
    ENDIF
    !
    IF (IERR .EQ. 0) THEN
       DO I=1,18
          NSDG(I)=2.3+REAL(I-1)*0.1
          NSFG(I)=3.3+REAL(I-1)*0.1
          DO J=1,18
             CS(I,J)=CS_DF(J,I)
             AL(I,J)=AL_DF(J,I)
          ENDDO
       ENDDO
       ! setup second derivative table for spline
       CALL SPLIE2(NSDG,NSFG,CS,18,18,CS2)
       CALL SPLIE2(NSDG,NSFG,AL,18,18,AL2)
       ! run bicubic spline interpolation
       CALL SPLIN2(NSDG,NSFG,CS,CS2,18,18,NSD,NSF,SIGMA)
       CALL SPLIN2(NSDG,NSFG,AL,AL2,18,18,NSD,NSF,ALPHA)
       !
    ENDIF
    !
  END SUBROUTINE TRANSITION_DF
  !---------------------------
  ! Subroutine RADIATIVE_DAMPING
  !---------------------------
  PURE SUBROUTINE RADIATIVE_DAMPING(LAMBDA,GAMMA_RAD)
    !
    IMPLICIT NONE
    !
    REAL(DP),   INTENT(IN)      :: LAMBDA
    REAL(SP),   INTENT(OUT)     :: GAMMA_RAD
    ! Calculate radiative damping assuming a classical oscillator
    ! LAMBDA must be in cm for GAMMA_RAD to be in s^(-1)
    GAMMA_RAD = 0.2223/LAMBDA**2.
    !
  END SUBROUTINE RADIATIVE_DAMPING
  !
  

  
  
!!$  !!**********************************************************************
!!$!! Computes the linewidth (half width at half maxima) given
!!$!! the cross-section and velocity parameter and tempertature
!!$!! Paul Barklem Dec 1997
!!$!!**********************************************************************
!!$!
!!$!! CROSS = cross section in atomi!! units at 10000m/s
!!$!! ALPHA = velocity parameter
!!$!! TEMP = temperature in Kelvin
!!$!! HALFWIDTH = half width per unit hydrogen atom density (rad /s cm^-3)
!!$!! A = Atomic Mass of the absorbing atom eg. A(Mg)=24.32
!!$!! RMASS = reduced mass
!!$!! MEANVEL = mean kinetic particle velocity in m/s
!!$!!
!!$
!!$   SUBROUTINE GETWIDTH(CROSS,ALPHA,TEMP,A,HALFWIDTH)
!!$   USE BASTYPE
!!$   IMPLICIT NONE
!!$   REAL(DP), INTENT(IN)  :: CROSS,ALPHA,TEMP
!!$   REAL(DP), INTENT(OUT) :: HALFWIDTH
!!$!! local
!!$   REAL(DP) :: MEANVEL,RMASS,CROSSMEAN,CROSSM,A
!!$!!
!!$   REAL(DP), parameter :: K =1.380658E-23      !boltzmanns constant J/K
!!$   REAL(DP), parameter :: M0=1.660540E-27      !unit atomic mass kg (Carbon 12)
!!$   REAL(DP), parameter :: A0=5.29177249E-11    !bohr radius m
!!$!! 
!!$   CROSSM=CROSS*A0*A0                    !cross section to m^2
!!$   RMASS=M0/(ONE/1.008_DP+ONE/A)         !calculate reduced mass
!!$   MEANVEL=SQRT(8.0_DP*K*TEMP/PI/RMASS)  !array of mean velocities
!!$!!
!!$!! work out cross section at mean perturber velocity for this Temp
!!$!!
!!$   CROSSMEAN= CROSSM*((MEANVEL/10000.0_DP)**(-ALPHA)) 
!!$!!
!!$!! the half-half width per unit perturber density m^3 rad/s
!!$!!
!!$   HALFWIDTH=(FOUR/PI)**(ALPHA*OO2)*(DGAMMA((FOUR-ALPHA)*OO2))*MEANVEL*CROSSMEAN
!!$!!
!!$   HALFWIDTH=HALFWIDTH * 1.E6_DP ! to rad/s cm^3
!!$
!!$   END SUBROUTINE GETWIDTH
!!$
!!$
!!$!**********************************************************************
!!$!  The next routine is from the NETLIB archive
!!$!   Computes the Gamma function
!!$!**********************************************************************
!!$
!!$      DOUBLE PRECISION FUNCTION DGAMMA(X)
!!$      USE BASTYPE
!!$      IMPLICIT NONE
!!$!!---------------------------------------------------------------------
!!$!!    From NETLIB, group SPECFUN
!!$!!     'AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL FUNCTIONS',
!!$!!     LECTURE NOTES IN MATHEMATICS, 506, NUMERICAL ANALYSIS DUNDEE,
!!$!!     1975, G. A. WATSON (ED.), SPRINGER VERLAG, BERLIN, 1976.
!!$!!******************************************************************
!!$!!
!!$!! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!!$!!
!!$!! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!!$!!          1.0 + EPS .GT. 1.0
!!$!! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!!$!!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!!$!!                  GAMMA(XBIG) = XINF.
!!$!! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!!$!!          1/XMININ IS MACHINE REPRESENTABLE.
!!$!! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER.
!!$!!
!!$!!
!!$!! ERROR RETURNS
!!$!!
!!$!!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!!$!!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!!$!!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!!$!!
!!$!! OTHER SUBPROGRAMS REQUIRED (DOUBLE PRECISION VERSION)
!!$!!
!!$!!     DBLE,DEXP,DLOG,DSIN,FLOAT,IFIX,SNGL
!!$!!
!!$      REAL(DP) :: EPS,FACT,HALF,RES,SQRTPI
!!$      REAL(DP) :: SUM,TWELVE,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z
!!$      INTEGER(I4B) :: I,J,N
!!$      LOGICAL(LGT) :: PARITY
!!$!!---------------------------------------------------------------------
!!$!!  MATHEMATICAL CONSTANTS
!!$!!---------------------------------------------------------------------
!!$      DATA HALF,TWELVE/0.5D0,12.0D0/
!!$      DATA SQRTPI/0.9189385332046727417803297D0/
!!$!!---------------------------------------------------------------------
!!$!!  MACHINE DEPENDENT PARAMETERS
!!$!!---------------------------------------------------------------------
!!$      DATA XBIG,XMININ,EPS/34.844D0,5.883D-39,2.776D-17/
!!$      DATA XINF/1.7014D38/
!!$!!---------------------------------------------------------------------
!!$!!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!!$!!     APPROXIMATION OVER (1,2).
!!$!!---------------------------------------------------------------------
!!$    REAL(DP), PARAMETER :: P(8) = &
!!$      (/-1.71618513886549492533811D+0, 2.47656508055759199108314D+1, &
!!$        -3.79804256470945635097577D+2, 6.29331155312818442661052D+2, &
!!$         8.66966202790413211295064D+2,-3.14512729688483675254357D+4, &
!!$        -3.61444134186911729807069D+4, 6.64561438202405440627855D+4/)
!!$    REAL(DP), PARAMETER :: Q(8) = &
!!$       (/-3.08402300119738975254353D+1, 3.15350626979604161529144D+2,&
!!$         -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,&
!!$          2.25381184209801510330112D+4, 4.75584627752788110767815D+3,&
!!$         -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/)
!!$!!---------------------------------------------------------------------
!!$!!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!!$!!---------------------------------------------------------------------
!!$     REAL(DP), PARAMETER :: C(8) = &
!!$       (/-1.910444077728D-03,   8.4171387781295D-04, &
!!$         -5.952379913043012D-04,7.93650793500350248D-04, &
!!$         -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
!!$          5.7083835261D-03, 0.D0/)
!!$!!---------------------------------------------------------------------
!!$      PARITY = .FALSE.
!!$      FACT = ONE
!!$      N = 0
!!$      Y = X
!!$      IF (Y .GT. ZERO) GO TO 200
!!$!!---------------------------------------------------------------------
!!$!!  ARGUMENT IS NEGATIVE
!!$!!---------------------------------------------------------------------
!!$      Y = -X
!!$      J = IFIX(SNGL(Y))
!!$      RES = Y - DBLE(FLOAT(J))
!!$      IF (RES .EQ. ZERO) GO TO 700
!!$      IF (J .NE. (J/2)*2) PARITY = .TRUE.
!!$      FACT = -PI / DSIN(PI*RES)
!!$      Y = Y + ONE
!!$!!---------------------------------------------------------------------
!!$!!  ARGUMENT IS POSITIVE
!!$!!---------------------------------------------------------------------
!!$  200 IF (Y .LT. EPS) GO TO 650
!!$      IF (Y .GE. TWELVE) GO TO 300
!!$      Y1 = Y
!!$      IF (Y .GE. ONE) GO TO 210
!!$!!---------------------------------------------------------------------
!!$!!  0.0 .LT. ARGUMENT .LT. 1.0
!!$!!---------------------------------------------------------------------
!!$      Z = Y
!!$      Y = Y + ONE
!!$      GO TO 250
!!$!!---------------------------------------------------------------------
!!$!!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!!$!!---------------------------------------------------------------------
!!$  210 N = IFIX(SNGL(Y)) - 1
!!$      Y = Y - DBLE(FLOAT(N))
!!$      Z = Y - ONE
!!$!!---------------------------------------------------------------------
!!$!!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!!$!!---------------------------------------------------------------------
!!$  250 XNUM = ZERO
!!$      XDEN = ONE
!!$      DO 260 I = 1, 8
!!$         XNUM = (XNUM + P(I)) * Z
!!$         XDEN = XDEN * Z + Q(I)
!!$  260 CONTINUE
!!$      RES = XNUM / XDEN + ONE
!!$      IF (Y .EQ. Y1) GO TO 900
!!$      IF (Y1 .GT. Y) GO TO 280
!!$!!---------------------------------------------------------------------
!!$!!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!!$!!---------------------------------------------------------------------
!!$      RES = RES / Y1
!!$      GO TO 900
!!$!!---------------------------------------------------------------------
!!$!!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!!$!!---------------------------------------------------------------------
!!$  280 DO 290 I = 1, N
!!$         RES = RES * Y
!!$         Y = Y + ONE
!!$  290 CONTINUE
!!$      GO TO 900
!!$!!---------------------------------------------------------------------
!!$!!  EVALUATE FOR ARGUMENT .GE. 12.0,
!!$!!---------------------------------------------------------------------
!!$  300 IF (Y .GT. XBIG) GO TO 700
!!$      YSQ = Y * Y
!!$      SUM = C(7)
!!$      DO 350 I = 1, 6
!!$         SUM = SUM / YSQ + C(I)
!!$  350 CONTINUE
!!$      SUM = SUM/Y - Y + SQRTPI
!!$      SUM = SUM + (Y-HALF)*DLOG(Y)
!!$      RES = DEXP(SUM)
!!$      GO TO 900
!!$!!---------------------------------------------------------------------
!!$!!  ARGUMENT .LT. EPS
!!$!!---------------------------------------------------------------------
!!$  650 IF (Y .LT. XMININ) GO TO 700
!!$      RES = ONE / Y
!!$      GO TO 900
!!$!!---------------------------------------------------------------------
!!$!!  RETURN FOR SINGULARITIES, EXTREME ARGUMENTS, ETC.
!!$!!---------------------------------------------------------------------
!!$  700 DGAMMA = XINF
!!$      WRITE (IUC,*) 'function DGAMMA has extreme argument :'
!!$      WRITE (IUC,*) 'argument is : ',X
!!$      WRITE (IUC,*) 'function returns machine upper bound '
!!$      GO TO 950
!!$!!---------------------------------------------------------------------
!!$!!  FINAL ADJUSTMENTS AND RETURN
!!$!!---------------------------------------------------------------------
!!$  900 IF (PARITY) RES = -RES
!!$      IF (FACT .NE. ONE) RES = FACT / RES
!!$      DGAMMA = RES
!!$  950 RETURN
!!$!! ---------- LAST CARD OF GAMMA ----------
!!$      END FUNCTION DGAMMA











END MODULE DAMPING
