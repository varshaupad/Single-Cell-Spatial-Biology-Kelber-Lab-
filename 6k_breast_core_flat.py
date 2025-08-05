# """
# Author: Varsha Upadhyayulla
# Project Title: 6K Breast Core Groupings
# Date Created: 8/16/24
# Date Last Modified: 12/20/24
# """
import pandas as pd

# core_names = {1:'F8', 2:'E8', 3:'D8', 4:'C8', 5:'G8', 6:'G7', 
#                 7:'F7', 8:'E7', 9:'D7', 10:'C7', 11:'C6', 12:'D6',
#                 13:'E6', 14:'F6', 15:'G6', 16:'G5', 17:'F5', 18:'E5',
#                 19:'D5', 20:'C5', 21:'G4', 22:'F4', 23:'E4', 24:'D4', 
#                 25:'C4', 26:'C3', 27:'D3', 28:'E3', 29:'F3', 30:'G3',
#                 31:'G2', 32:'F2', 33:'E2', 34:'D2', 35:'C2', 36:'C1', 
#                 37:'D1', 38:'E1', 39:'F1', 40:'G1', 41:'C9', 42:'D9', 
#                 43:'E9', 44:'F9', 45:'G9', 46:'H9', 47:'H8', 48:'H7',
#                 49:'H6', 50:'H5', 51:'H4', 52:'H3', 53:'H2', 54:'H1'}

# range_to_core_sample = {(1, 9): 1, (10, 18): 2, (19, 27): 3, (28, 36): 4, (37, 45): 5, 
#                             (46, 54): 6, (55, 63): 7, (64, 72): 8, (73, 81): 9, (82, 90): 10,
#                             (91, 99): 11, (100, 108): 12, (109, 117): 13, (118, 126): 14, (127, 135): 15,
#                             (136, 144): 16, (145, 153): 17, (154, 162): 18, (163, 171): 19, (172, 180): 20,
#                             (181, 189): 21, (190, 198): 22, (199, 207): 23, (208, 216): 24, (217, 225): 25, 
#                             (226, 234): 26, (235, 243): 27, (244, 252): 28, (253, 261): 29, (262, 270): 30,
#                             (271, 279): 31, (280, 288): 32, (289, 297): 33, (298, 306): 34, (307, 315): 35,
#                             (316, 324): 36, (325, 333): 37, (334, 342): 38, (343, 351): 39, (352, 360): 40,
#                             (361, 362): 41, (363, 364): 42, (365, 366): 43, (367, 368): 44, (369, 370): 45,
#                             (371, 371): 46, (372, 373): 47, (374, 375): 48, (376, 377): 49, (378, 379): 50,
#                             (380, 381): 51, (382, 383): 52, (384, 385): 53, (386, 387): 54}

# groupings = {'HER2': ('C1', 'C6', 'C8', 'C9', 'D1', 'D6', 'D8', 'E5', 
#                       'E7', 'F5', 'F7', 'G1', 'G2', 'G3', 'G6', 'G8', 
#                       'H1', 'H2', 'H3', 'H6', 'H8'),

#              'TNBC': ('C7', 'D7', 'E1', 'E2', 'E3', 'F1', 'G5', 'G7', 'G9', 
#                      'H5', 'H7', 'H9'),    

#              'ER/PR': ('C3','C4', 'C5', 'D2', 'D3', 'D4', 'D5', 'E4', 'E6', 'E8', 'E9', 
#                        'F2', 'F3', 'F4', 'F6', 'F8', 'F9', 'G4', 'H4'),    

#              'Primary': ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 
#                          'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 
#                          'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9',),
                         
#              'Met': ('D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 
#                      'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 
#                      'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9'),

#              'Over or Equal to 40': ('C1', 'C2', 'C3', 'C5', 'C7', 'C8', 'C9',
#                                      'D1', 'D2', 'D3', 'D5', 'D7', 'D8', 'D9',
#                                      'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8',
#                                      'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8',
#                                      'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G8', 'G9',
#                                      'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H8', 'H9'),

#              'Under 40': ('C4', 'C6', 'D4', 'D6', 'E9', 'F9', 'G7', 'H7')}

# her2_primary_cores = []
# her2_met_cores = []
# tnbc_primary_cores = []
# tnbc_met_cores = []
# erpr_primary_cores = []
# erpr_met_cores = []
# over_40_primary_cores = []
# over_40_met_cores = []
# under_40_primary_cores = []
# under_40_met_cores = []


with open("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Primary.csv", "rt") as file1, \
    open("/data/kelberj_shared/data/Brm961b_022-6KPanel/grouped/Breast_Core_Groups/6K_RNA_HER2_Met.csv", "rt") as file2:
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    median1 = df1['PTPRC'].median()
    median2 = df2['PTPRC'].median()

    low_ptprc_cells1 = df1[df1['PTPRC'] <= median1]
    high_ptprc_cells1 = df1[df1['PTPRC'] > median1]

    low_ptprc_cells2 = df2[df2['PTPRC'] <= median2]
    high_ptprpc_cells2 = df2[df2['PTPRC'] > median2]

    with open('6K_RNA_HER2_Primary_Low_PTPRC.csv', mode='w') as file3, \
         open('6K_RNA_HER2_Primary_High_PTPRC.csv', mode='w') as file4, \
         open('6K_RNA_HER2_Met_Low_PTPRC.csv', mode='w') as file5, \
         open('6K_RNA_HER2_Met_High_PTPRC.csv', mode='w') as file6:

        low_ptprc_cells1.to_csv(file3, index = False)
        high_ptprc_cells1.to_csv(file4, index = False)

        low_ptprc_cells2.to_csv(file5, index = False)
        high_ptprpc_cells2.to_csv(file6, index = False)
        
# with open("C:\\Users\\Raynah_Cheng1\\OneDrive - Baylor University\\Spatial_Biology\\All_Data\\6K_RNA_Data\\BRM961B_022_6KPanel_exprMat_file.csv", "r") as file:

     full_df = pd.read_csv(file)

     def map_to_core_sample(fov):
         for (start, end), core_sample in range_to_core_sample.items():
             if start <= fov <= end:
                 return core_sample
           
     full_df['core sample'] = full_df['fov'].apply(map_to_core_sample)
     full_df['core sample'] = full_df['core sample'].map(core_names)
    
     column_titles = full_df.columns.tolist()

     grouped = full_df.groupby('core sample')

#     with open('6K RNA HER2 Primary.csv', mode='w') as file1, \
#         open('6K RNA HER2 Met.csv', mode='w') as file2, \
#         open('6K RNA TNBC Primary.csv', mode='w') as file3, \
#         open('6K RNA TNBC Met.csv', mode='w') as file4, \
#         open('6K RNA ERPR Primary.csv', mode='w') as file5, \
#         open('6K RNA ERPR Met.csv', mode='w') as file6, \
#         open('6K RNA Over or Equal to 40 Primary.csv', mode='w') as file7, \
#         open('6K RNA Over or Equal to 40 Met.csv', mode='w') as file8, \
#         open('6K RNA Under 40 Primary.csv', mode='w') as file9, \
#         open('6K RNA Under 40 Met.csv', mode='w') as file10:

#         # Write the headers to each file
#         file1.write(','.join(column_titles) + '\n')
#         file2.write(','.join(column_titles) + '\n')
#         file3.write(','.join(column_titles) + '\n')
#         file4.write(','.join(column_titles) + '\n')
#         file5.write(','.join(column_titles) + '\n')
#         file6.write(','.join(column_titles) + '\n')
#         file7.write(','.join(column_titles) + '\n')
#         file8.write(','.join(column_titles) + '\n')
#         file9.write(','.join(column_titles) + '\n')
#         file10.write(','.join(column_titles) + '\n')

     for core_sample, group in grouped:
         core_name = core_sample

         if core_name in groupings['HER2'] and core_name in groupings['Primary']:
             her2_primary_cores.append(core_name)
             group.to_csv('6K RNA HER2 Primary.csv', mode = "a", index = False, header = False)
        
         if core_name in groupings['HER2'] and core_name in groupings['Met']:
            
             her2_met_cores.append(core_name)
             group.to_csv('6K RNA HER2 Met.csv', mode = "a", index = False, header = False)

         if core_name in groupings['TNBC'] and core_name in groupings['Primary']:
            
             tnbc_primary_cores.append(core_name)
             group.to_csv('6K RNA TNBC Primary.csv', mode = "a", index = False, header = False)

         if core_name in groupings['TNBC'] and core_name in groupings['Met']:
            
             tnbc_met_cores.append(core_name)
             group.to_csv('6K RNA TNBC Met.csv', mode = "a", index = False, header = False)

         if core_name in groupings['ER/PR'] and core_name in groupings['Primary']:
            
             erpr_primary_cores.append(core_name)
             group.to_csv('6K RNA ERPR Primary.csv', mode = "a", index = False, header = False)

         if core_name in groupings['ER/PR'] and core_name in groupings['Met']:
            
             erpr_met_cores.append(core_name)
             group.to_csv('6K RNA ERPR Met.csv', mode = "a", index = False, header = False)

#         if core_name in groupings['Over or Equal to 40'] and core_name in groupings['Primary']:
            
#             over_40_primary_cores.append(core_name)
#             group.to_csv('6K RNA Over or Equal to 40 Primary.csv', mode = 'a', index = False, header = False)

#         if core_name in groupings['Over or Equal to 40'] and core_name in groupings['Met']:
            
#             over_40_met_cores.append(core_name)
#             group.to_csv('6K RNA Over or Equal to 40 Met.csv', mode = "a", index = False, header = False)

#         if core_name in groupings['Under 40'] and core_name in groupings['Primary']:
            
#             under_40_primary_cores.append(core_name)
#             group.to_csv('6K RNA Under 40 Primary.csv', mode = "a", index = False, header = False)

#         if core_name in groupings['Under 40'] and core_name in groupings['Met']:
            
#             under_40_met_cores.append(core_name)
#             group.to_csv('6K RNA Under 40 Met.csv', mode = "a", index = False, header = False)
