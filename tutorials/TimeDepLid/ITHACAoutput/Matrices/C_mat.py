C=np.zeros([12,12,12])
C[0,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0]])

C[1,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-1.139225963e-08,-0.001025271579,0.001108249418,-0.0004468746192,-0.0003078971081,0.001207283398,-0.001190571406,0.001937176906,-0.0002454491077,0.001210109827,0.0004960098441],
[0,-5.278833364e-06,0.3508712157,-0.4548582632,-0.7504325124,0.6453149977,-0.3968979671,-0.4973497426,0.1657734961,-0.4791988475,0.1338193445,-0.3745293326],
[0,2.43919222e-06,-0.3142052338,0.6504734007,-0.8040636972,-0.7046357594,0.8407216151,-0.09641546306,0.4590964238,-0.03207412633,0.4867727961,0.1732772146],
[0,-5.714187879e-06,0.3075625128,-0.4181330962,0.5243387053,-0.2787246301,-1.592510042,0.880422546,-0.1000093068,-0.5659750196,-0.01057022378,-0.3337677484],
[0,-4.885363986e-06,-0.003967138176,0.1863667666,-0.6191237343,0.8281006727,-0.4428271644,-2.045361078,1.563104496,-0.2039978216,0.3592768366,0.1151028666],
[0,7.531215484e-06,0.220443093,0.3760797641,-0.2548304319,-1.010510085,0.8283393726,0.2827998619,-2.046730833,1.315666673,-0.2210410005,-0.4213624617],
[0,-1.490580649e-05,-0.1735585984,0.04153490408,0.7546598164,-0.4051777165,-0.3850071734,0.9260578431,-0.6953507365,-1.480956091,1.299332513,-0.3847366457],
[0,1.688660217e-05,0.3553489504,0.04906571021,-0.1668914284,0.6126747668,-0.6034773044,-1.129501526,1.039918453,-0.5915724331,-1.188360566,1.385214199],
[0,-9.013996318e-06,0.006712576553,0.1966437014,-0.2899339602,-0.2818773136,1.170687705,-0.3071025809,-1.090516471,1.500363139,-1.427312668,-0.9753388156],
[0,1.085980257e-05,0.1696512543,0.1025636232,0.2147125785,-0.006593160926,-0.2909788429,0.8122144892,-0.2040425253,-1.674445156,1.84028621,-1.218406502],
[0,-2.371790558e-06,0.09622456222,0.108876819,-0.1875917672,0.1337370193,-0.1112864359,-0.5224597561,1.143352448,-0.04472651618,-1.602461944,1.878159997]])

C[2,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0.001017941709,-0.001715248545,2.415208934,-0.4418181735,-0.7314852085,0.8689888831,-0.4639480646,0.7992354663,0.1058584063,0.3869521512,0.3321769698],
[0,-0.35341527,-0.6054216465,-558.3210484,-29.96293741,401.8744393,12.77417616,-141.8396076,220.7446826,-119.0140071,60.77089571,-2.941710023],
[0,0.3161444962,-0.4644412059,700.3511576,-227.3974416,76.84189464,723.928766,-148.6473358,275.1120206,355.7863789,128.827727,227.067865],
[0,-0.3140747107,-0.6372334659,-737.6400939,-111.2878005,342.0253849,-14.02566299,363.0799635,-31.9085477,-425.4230027,220.0197649,-61.69296484],
[0,0.00109156098,-0.3038162249,676.5798854,-244.4747351,-455.7557527,239.9764875,-35.87711296,825.611166,-87.95565511,-39.54283613,344.890349],
[0,-0.2154332871,0.576666822,-633.4222696,276.9959123,240.2130776,-397.987041,169.1315486,74.05055154,740.3177722,-403.989819,-319.2195004],
[0,0.1613912945,-2.116998795,69.9156257,-77.02213341,438.3947691,-61.03758576,-531.504024,65.40249015,140.9596786,796.5070849,-406.3582575],
[0,-0.3431699084,1.959326622,-382.2422576,-191.0715562,16.59946859,594.5236931,115.767603,-494.1037307,231.2413647,195.5764934,699.2149658],
[0,-0.01316309697,-1.229195603,-107.5650978,-82.7770545,-84.88933126,-3.702296733,509.8202449,64.20788692,-405.0702186,-173.7711454,483.5827813],
[0,-0.1643304931,0.7284340496,-181.0971259,48.47759138,78.32280473,-261.409347,-71.14072805,621.0996939,174.176317,-533.0272571,-54.29456624],
[0,-0.09676506783,0.1054561617,-141.3342309,-108.6911801,106.5742942,185.3836467,-303.6522086,-133.2360884,824.1677831,-45.05203891,-445.4875575]])

C[3,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-0.001118239583,-2.417428414,-0.00268181664,3.962316265,-1.531985822,0.3863554136,-0.07573962783,0.2865622871,0.2855558949,0.06074516869,0.2329141842],
[0,0.4487843841,557.7553174,-2.047375908,-2291.023256,156.7808363,201.8541824,-167.370093,139.6300902,41.01060858,-221.667142,-28.03089675],
[0,-0.648718491,-699.1534668,0.09930876178,667.6024921,-1519.59944,201.1091688,722.1661338,-204.2472807,69.47978283,175.2191954,-127.7818315],
[0,0.41729299,734.6805194,1.178564443,-754.6696343,520.1796217,-1180.010867,108.3947735,905.859382,-239.0800668,-77.22160526,342.1194541],
[0,-0.195736174,-677.9620737,-3.710058357,323.1010268,-82.21423337,104.246387,-1365.831468,502.6440146,1068.629803,-467.4195146,136.2200122],
[0,-0.366561399,636.7552409,2.535543765,-937.7534513,-872.9045385,380.4651301,357.8204343,-1837.955233,612.3050259,1165.85026,-979.6914988],
[0,-0.04841371922,-73.97002751,0.205385546,700.9453122,-58.7065471,-396.2214034,467.8020784,-267.1572174,-1572.388873,968.0257579,835.8027269],
[0,-0.04184998505,385.2957113,0.03964233757,-645.8927854,359.4037623,-129.3141455,-966.2059995,546.5728201,25.84729953,-1526.995075,1210.044394],
[0,-0.1962761013,107.381289,1.046879245,-522.1463617,-284.1054925,1111.37294,25.19281375,-1108.409505,955.9062646,-602.3446252,-1754.83138],
[0,-0.09960033843,180.9930965,0.09496217698,-139.9977894,-209.7239566,-99.96592311,853.8933423,-3.846073244,-1338.643878,1165.962842,-371.5516005],
[0,-0.1074040384,142.2807651,0.2850813999,-362.6442643,-57.56684892,-162.457475,-284.1215865,960.0040428,251.9314561,-1465.931082,1239.963079]])

C[4,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0.0004452834482,0.4416693265,-3.96150712,-9.286488341e-05,8.224374827,0.06013109654,-3.819047111,0.8330275743,-0.6487020249,0.01590021095,-0.1662138724],
[0,0.7518619516,29.63120898,2292.915307,-0.7621715064,-2860.761348,-1293.890929,1159.263003,686.7464327,397.5432197,488.7115938,168.7557072],
[0,0.7990765511,227.1952562,-670.038799,1.40349768,2479.340516,-1147.264851,-2055.944308,456.907806,-21.87371227,-37.02329045,89.80331108],
[0,-0.5243471829,111.8214652,753.8438344,-0.7485016943,-1975.228192,382.8393027,576.5651733,-1100.904138,201.1329733,289.4069957,-513.7027772],
[0,0.6255373393,244.6757528,-317.097721,-2.316390934,2031.501776,46.76940215,-297.4779112,353.6359248,-479.3723974,934.3839617,264.4821483],
[0,0.2451462607,-278.6964024,932.0417488,2.668776601,-1667.04073,-270.6325013,740.8345143,-81.65206604,-863.4486809,-127.4637957,1357.322051],
[0,-0.7590601397,77.27771807,-704.2867592,0.6544510577,-65.52712633,-162.7962908,6.332383765,436.8769318,-211.7375814,-1107.837879,-196.1398271],
[0,0.1737146712,192.5214381,650.5816147,-1.594412574,-913.082032,-464.5453261,261.4667537,-140.6852652,530.7027237,-107.892078,-1137.854657],
[0,0.2782605934,79.21140376,517.4675577,1.792631288,-196.1778458,-322.5545602,-42.12331532,102.5320184,-416.3655032,1112.589936,-47.02946589],
[0,-0.2083337919,-45.75368761,141.0797603,0.5433540111,-628.7136656,185.3430667,37.20260025,-291.4883839,221.7490004,-401.2212993,981.6472909],
[0,0.1830237843,107.0620257,362.7276199,-1.117584749,-466.0857635,-199.7831675,668.4819759,-217.8833206,-521.7512707,385.7431103,-439.8757414]])

C[5,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0.0002958656195,0.7266031353,1.526680441,-8.224992609,-0.003540338979,8.328352357,0.2986208229,0.1050765969,1.564135122,0.01202393219,1.153283307],
[0,-0.6494441567,-403.7661485,-158.9974824,2861.175292,-1.353413819,-2329.754387,-2280.855555,100.4949377,273.0508996,72.38155126,37.74031073],
[0,0.7110530698,-76.22431544,1522.578727,-2482.219539,2.705594446,3781.918933,-613.3808354,-714.452198,982.1635081,201.0062786,863.267671],
[0,0.2654148188,-345.3882505,-525.5564622,1977.9125,-4.936805619,-2673.825589,529.8398811,224.428047,-2017.67544,45.00762982,-247.686587],
[0,-0.8321176408,454.1528512,77.60681649,-2029.234249,-0.941723393,1539.25386,-1320.964239,1020.893967,469.2767337,-923.9547177,628.1129782],
[0,1.026277892,-236.6046483,883.0929127,1659.258054,6.583113376,-782.7918022,336.0487732,-351.0424019,531.4994967,-589.3134607,-410.3963643],
[0,0.3808171409,-446.3095916,50.37080454,66.3234712,-9.23499656,-312.8375472,337.6226836,343.9085435,-323.9286676,918.5141811,-888.0476701],
[0,-0.5908754058,-8.618720441,-354.8974178,915.7122013,7.911016159,-631.1776427,-788.4676691,143.7220936,135.0738536,-315.0758599,659.1276751],
[0,0.2800408236,82.88040413,289.482682,186.226846,1.566079654,-44.03119574,-118.5754534,-308.7615745,455.0119369,-509.610795,460.0592964],
[0,0.01084418604,-77.77477786,207.3802549,635.126064,-1.888696971,-655.9237628,-238.3636346,402.9801675,-505.2450222,405.4126511,-285.2242909],
[0,-0.1277672336,-103.9410472,60.99544477,462.9855596,5.518975267,-51.64382489,-684.4227372,177.3644443,624.7558599,-799.8817849,645.6240642]])

C[6,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-0.001204314308,-0.8659936466,-0.3855052182,-0.05855598773,-8.325999465,4.178461361e-05,14.77767081,-0.9574243237,-5.820751696,1.512408767,-1.963975161],
[0,0.3959959378,-12.0143387,-202.4617212,1295.133852,2329.352357,0.149900788,-4693.851936,-2918.7962,1477.704876,763.8299993,21.12060447],
[0,-0.8448492756,-724.483099,-202.9456226,1148.007681,-3784.97848,1.551031782,5866.575375,-1438.837535,-3583.05085,743.776255,-302.5273224],
[0,1.602519868,15.84739693,1185.039912,-384.553976,2679.896795,-1.492086537,-4213.083744,2234.943925,2739.753056,-2426.124608,795.7388913],
[0,0.4342774192,-239.9044435,-108.2244916,-42.81749194,-1544.167049,1.401174712,2482.849142,-2345.6291,56.34639771,859.9479125,-2659.395524],
[0,-0.8303447432,396.828323,-383.3355207,271.4881318,781.3223897,-0.753280999,-2112.512878,-436.2026761,-185.8182506,561.2440498,-423.5455778],
[0,0.4052054599,66.87856638,405.5348369,158.571333,327.6805033,-2.943154146,384.4024788,947.6155158,295.8844233,-211.2661453,1482.661599],
[0,0.5786136466,-600.6641763,119.6040565,467.404005,614.7448408,6.086889821,-1490.806151,-480.0366212,824.8128974,-394.6385212,-229.0738813],
[0,-1.162672182,6.731059247,-1111.845057,325.3690858,47.76104505,-2.96465761,-650.4882488,-739.1094168,253.273849,224.7457642,-1212.816362],
[0,0.2880406864,259.7061855,103.3588872,-190.711139,655.0986047,2.152738322,-884.301489,-234.6485538,337.7291775,-12.99940923,343.3428943],
[0,0.1014355096,-188.1319597,155.526514,206.1002803,41.51610102,1.873328343,-704.3839237,-292.9779299,436.1306831,-78.04750165,-180.0244091]])

C[7,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0.001178968678,0.4602193429,0.07273603881,3.817903214,-0.3049693776,-14.77397751,-0.002862820293,12.32394907,-1.928662407,-0.6156890451,1.166253574],
[0,0.4944575699,140.4454729,167.0699816,-1160.3213,2279.787827,4694.674264,-1.419765318,-2925.707141,-3001.833643,1009.472578,676.4454014],
[0,0.09260295346,148.5394054,-724.1669372,2057.385455,611.3169102,-5866.807216,0.7685449054,5280.133914,-1837.64202,-1718.673204,283.856384],
[0,-0.8824544875,-365.2031951,-107.7406059,-580.0367875,-530.6205411,4214.814455,-3.170962651,-3758.854211,2616.193422,1376.424938,-2530.910067],
[0,2.048414109,35.33201681,1368.852246,295.7897041,1322.305366,-2482.148132,-0.3455962537,2706.038668,-2360.126425,2591.064184,1279.852157],
[0,-0.2926936603,-170.475611,-363.3421061,-735.7859153,-343.0004672,2112.852305,1.430986923,-2193.255327,-142.8991924,-1158.232089,255.4204835],
[0,-0.9404985332,526.1512894,-472.8808602,-8.589807256,-344.3359936,-381.803544,-5.561448217,-227.5369711,144.7839945,109.7015766,-977.9822589],
[0,1.153129519,-108.0891793,974.952696,-261.4075806,803.0670106,1484.611121,5.128694885,-646.0541958,-122.6494576,493.917506,251.537752],
[0,0.2749037139,-517.982774,-36.5856941,43.60710101,100.6907329,660.8876763,-5.648311682,-652.9030566,-258.934842,418.1642169,-112.2082372],
[0,-0.7866009452,78.64268238,-846.9594154,-36.68717023,252.537072,875.2326617,6.826094044,-1041.685526,-337.4848957,84.10909606,-56.13959272],
[0,0.5118809854,299.8798429,284.6997923,-672.2687327,679.3160828,709.302661,-4.391181727,-413.7889318,-543.7696577,404.1413857,-123.5158043]])

C[8,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-0.001929393553,-0.7963533942,-0.2853600958,-0.8312168781,-0.1011890501,0.9540315358,-12.31858433,-0.002083478236,20.4107698,-5.11660744,-3.078505783],
[0,-0.1636301932,-219.6684011,-139.880217,-685.3885994,-99.75526101,2917.977808,2928.42039,-1.385081961,-4843.185398,-3181.136734,2243.295128],
[0,-0.4547523158,-274.1101558,206.0255845,-458.1314124,717.8364233,1436.8363,-5280.957598,1.204903899,8561.639264,-2633.723289,-2630.878835],
[0,0.0992626871,32.67655266,-907.3406811,1103.757894,-226.982895,-2233.503256,3764.04338,-3.68998328,-6908.380173,3522.026333,2128.363678],
[0,-1.563018687,-825.615988,-503.8660096,-352.5687337,-1021.186661,2343.999305,-2705.2101,-0.2483205707,5169.677332,-4724.010104,2558.036784],
[0,2.060305677,-71.76940485,1845.854679,76.87933875,360.4246382,433.5401933,2190.69695,2.529767054,-2725.783976,1907.08872,-1031.82019],
[0,0.6981271594,-63.32634311,266.5847781,-431.9174983,-346.4016885,-945.4498494,237.2808594,-6.384851231,-1246.449559,907.8467747,386.3546707],
[0,-1.05182607,489.3422351,-551.0095685,139.1541493,-149.0234552,479.5309354,637.8422379,4.675415176,-1194.940813,-1314.892365,981.6909017],
[0,1.121488539,-54.69609739,1120.018206,-104.3779118,329.8247681,729.6575514,660.7822386,-1.35491955,-134.5512809,561.9136709,31.39356837],
[0,0.1730284432,-629.7443626,-6.970340205,292.6077186,-421.5323547,246.8032512,1030.082837,2.495343129,-1829.71066,369.6935137,931.1382555],
[0,-1.117229261,139.3729219,-952.2254716,216.9371242,-163.2629652,283.7670035,425.0329916,-4.200997009,-698.3915925,-1302.307848,1018.398822]])

C[9,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,0.0002321768081,-0.1082399756,-0.2892666033,0.6489160113,-1.570371863,5.825747773,1.922486011,-20.4061504,-0.0004906093099,19.29286519,-7.519461964],
[0,0.4718576224,117.3928296,-43.41806301,-397.3536257,-276.6958605,-1475.573652,2998.833127,4845.792157,-0.9605414824,-4747.118089,-3432.138586],
[0,0.02519774016,-356.8553364,-72.19601857,22.50908498,-985.9359853,3586.033819,1835.389845,-8561.765616,0.9912886158,8705.708878,-3624.584859],
[0,0.5673342629,425.0643673,241.3272994,-203.3144766,2019.416631,-2740.546203,-2618.507197,6913.406447,-2.840858512,-6965.875243,5690.930615],
[0,0.1904969883,86.1162192,-1074.309196,482.950492,-477.4808947,-51.83388811,2359.154053,-5170.674649,1.132844456,3446.361499,-6496.45158],
[0,-1.314000406,-740.6637944,-612.4827517,861.9264409,-529.5063655,183.5885679,140.8027132,2726.859129,-1.482911005,-1880.34029,1543.479502],
[0,1.484629582,-141.6235944,1578.280087,205.2664227,328.9379302,-296.188783,-153.2493421,1258.578747,-4.167388531,-393.3229419,2530.638434],
[0,0.5856900751,-231.6450648,-30.66061718,-524.3491363,-144.3770616,-820.0991165,131.5015608,1181.085351,6.263962048,-2234.671258,-1293.440324],
[0,-1.513466668,400.0678209,-960.4259274,414.2310761,-460.535976,-254.0200707,249.7828093,144.1566186,-7.765987776,-547.7031214,-765.4133146],
[0,1.690538296,-167.6249991,1344.706644,-219.8618517,515.7375286,-341.1976549,344.0710685,1825.616911,6.039830139,-1038.971687,1152.930638],
[0,0.01518887121,-831.9860658,-263.0850445,524.0518005,-642.0166963,-426.3231977,533.8188629,701.2698751,-2.126812645,-1352.341286,-648.6872888]])

C[10,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-0.001200740874,-0.3844207862,-0.05933752428,-0.01447349304,-0.007688130924,-1.516499029,0.622528466,5.109689686,-19.2907478,-0.0007942021048,27.03961982],
[0,-0.1273552513,-59.03756503,223.1804632,-488.0825535,-69.21706593,-766.3924401,-1005.543067,3177.629918,4749.049995,-0.8114557545,-5122.296792],
[0,-0.4815393734,-127.0162358,-174.0814293,37.04295602,-197.2516754,-746.870243,1721.349987,2632.012879,-8706.988009,0.9404037772,11766.59677],
[0,0.01184347676,-219.9893835,77.13056919,-289.0818606,-45.79678784,2427.70337,-1375.092832,-3525.241361,6971.378802,-3.42069213,-9692.179683],
[0,-0.349422596,41.58165336,470.5563595,-934.5290697,929.8138324,-865.0848277,-2586.761062,4722.106065,-3449.110272,2.022666055,7226.863057],
[0,0.2167001878,404.4231438,-1168.0483,129.6979369,587.2467633,-559.7823074,1158.749316,-1908.055392,1882.997518,-1.699146831,-4554.211768],
[0,-1.298972687,-796.1700455,-970.8056377,1110.592339,-921.1973964,212.7382002,-105.3506186,-917.198567,403.3192823,-6.041997486,-1321.449094],
[0,1.196770808,-193.9429683,1533.207975,102.4709499,324.3611139,389.0519206,-498.0418801,1325.437973,2222.063873,7.636891073,-657.7328964],
[0,1.43445221,176.7979146,603.7491475,-1109.408466,509.3302686,-222.8365863,-408.9078327,-572.7528664,562.2954198,-8.089255544,-1919.641072],
[0,-1.850724811,527.9502062,-1169.557495,398.8808883,-411.2971378,12.90456523,-91.04425275,-363.736637,1026.634644,6.261103042,-1929.892308],
[0,1.626663773,54.50679971,1474.438997,-382.7158537,816.6587225,67.59111417,-392.8851349,1298.357213,1358.580855,-2.544895213,-244.0676631]])

C[11,:,:]=np.array([[0,0,0,0,0,0,0,0,0,0,0,0],
[0,-0.0005100991088,-0.3347409916,-0.2368156834,0.1660316247,-1.159700368,1.968822655,-1.17348942,3.086163195,7.516162185,-27.03865695,0.001631898563],
[0,0.3651572374,0.6958989093,24.96025109,-168.3458023,-42.71906466,-18.23042322,-680.5295811,-2239.708622,3429.894748,5123.620879,-0.05054417333],
[0,-0.1779093893,-227.8527507,126.66418,-90.81679385,-864.9224329,304.2073492,-287.6269062,2634.76322,3623.116076,-11767.04469,1.18264162],
[0,0.3298233839,60.49469016,-342.7461441,513.294623,244.6897202,-793.2773509,2529.959349,-2128.101446,-5691.635583,9695.223784,-2.128424096],
[0,-0.1286968371,-347.8389171,-141.5259478,-262.6094263,-635.9097601,2663.649835,-1284.076953,-2554.542374,6495.492123,-7229.061291,1.93173706],
[0,0.4340822467,321.3915819,985.6917558,-1361.377949,420.3360042,418.2701943,-256.8084467,1036.034388,-1547.853637,4558.765677,-0.9623861953],
[0,0.3730689917,404.1026196,-838.1703991,195.9128845,881.9860237,-1476.718127,972.8570307,-382.6018416,-2531.846868,1325.818394,-1.698997915],
[0,-1.388676145,-699.5119276,-1215.000484,1141.742628,-662.9934491,228.3812752,-247.1950494,-988.6589618,1296.913973,651.01626,1.037429333],
[0,0.9821599397,-482.4026854,1761.105408,39.37301675,-449.7255091,1206.909361,103.0373428,-15.08138137,750.5427219,1934.548594,-3.009233235],
[0,1.213923809,53.57634111,368.700715,-978.8036068,276.4921678,-334.7398815,61.38483537,-942.9643707,-1136.236987,1913.234419,4.998764556],
[0,-1.886735376,440.8628722,-1242.602841,436.7976226,-649.8726109,178.3690086,116.3397531,-1010.724415,633.4202813,258.5249999,-7.603652067]])

