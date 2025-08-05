"""
Author: Raynah Cheng and Varsha Upadhyayulla
Assignment Title: 6K Pancreas Core Groupings
Date Created: 9/4/24
Date Last Modified: 12/20/24
"""
import pandas as pd

core_names = {1:'G8', 2:'F8', 3:'E8', 4:'D8', 5:'C8', 6:'B8', 
                7:'G7', 8:'F7', 9:'E7', 10:'D7', 11:'C7', 12:'B7',
                13:'G6', 14:'F6', 15:'E6', 16:'D6', 17:'C6', 18:'B6',
                19:'G5', 20:'F5', 21:'E5', 22:'D5', 23:'C5', 24:'B5', 
                25:'B4', 26:'C4', 27:'D4', 28:'E4', 29:'F4', 30:'G4',
                31:'G3', 32:'F3', 33:'E3', 34:'D3', 35:'C3', 36:'B3', 
                37:'G2', 38:'F2', 39:'E2', 40:'D2', 41:'C2', 42:'B2', 
                43:'G1', 44:'F1', 45:'E1', 46:'D1', 47:'C1', 48:'B1',
                49:'H1', 50:'H2', 51:'H3', 52:'H4', 53:'H5', 54:'H6',
                55: 'H7', 56: 'H8'}
    
range_to_core_sample = {(1, 6): 1, (7, 12): 2, (13, 18): 3, (19, 24): 4, (25, 30): 5,
                        (31, 36): 6, (37, 42): 7, (43, 48): 8, (49, 54): 9, (55, 60): 10,
                        (61, 66): 11, (67, 72): 12, (73, 78): 13, (79, 84): 14, (85, 90): 15, 
                        (91, 96): 16, (97, 102): 17, (103, 108): 18, (109, 114): 19, (115, 120): 20,
                        (121, 126): 21, (127, 132): 22, (133, 138): 23, (139, 144): 24, (145, 150): 25,
                        (151, 156): 26, (157, 162): 27, (163, 168): 28, (169, 174): 29, (175, 180): 30,
                        (181, 186): 31, (187, 192): 32, (193, 198): 33, (199, 204): 34, (205, 210): 35,
                        (211, 216): 36, (217, 222): 37, (223, 228): 38, (229, 234): 39, (235, 240): 40,
                        (241, 246): 41, (247, 252): 42, (253, 258): 43, (259, 264): 44, (265, 270): 45,
                        (271, 276): 46, (277, 282): 47, (283, 288): 48, (289, 290): 49, (291, 292): 50, 
                        (293, 294): 51, (295, 296): 52, (297, 298): 53, (299, 300): 54, (301, 302): 55,
                        (303, 304): 56}

groupings = {"Normal": ('H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8'),
             "NAT": ('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8'),
             "Grade 2": ('B1', 'B3', 'B4', 'B5', 'B6', 'B8', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'D1', 'D2'),
             "Grade 3": ('B7', 'D4', 'D5', 'D6', 'E1', 'E2', 'E4',  'E6', 'E7', 'E8', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8'),
             "Basal": ('E4', 'E5', 'E6', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'B7', 'B8', 'C1', 'C2', 'C3', 'C7', 'C8',
                       'D4', 'D5', 'D6', 'D7', 'D8'),
             "Classical": ('E1', 'E2', 'E3', 'E7', 'E8', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'C4', 'C5', 'C6', 'D1', 'D2', 'D3'),
             "Over or Equal to 60": ('B1', 'B2', 'B3', 'B7', 'B8', 'C1', 'C2', 'C3', 'C7', 'C8', 'D7', 'D8', 
                                     'E4', 'E5', 'E6', 'E7', 'E8', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'G4', 'G5', 'G6', 
                                     'G7', 'G8', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8'),
             "Under 60": ('B4', 'B5', 'B6', 'C4', 'C5', 'C6', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'E1', 'E2', 'E3', 'F7', 'F8',
                          'G1', 'G2', 'G3')}

with open("/data/kelberj_shared/data/PA722_151_6KPanel/flatfile/PA722_151_6KPanel_exprMat_file.csv", "r") as file:
    full_df = pd.read_csv(file)

    def map_to_core_sample(fov):
        for (start, end), core_sample in range_to_core_sample.items():
            if start <= fov <= end:
                return core_sample
            
    full_df['core sample'] = full_df['fov'].apply(map_to_core_sample)
    full_df['core sample'] = full_df['core sample'].map(core_names)

    column_titles = full_df.columns.tolist()


    grouped = full_df.groupby('core sample')

    with open('6K_RNA_Over60_Malignant.csv', mode='w') as file1, \
        open('6K_RNA_Under60_Malignant.csv', mode='w') as file2:

        # Write the headers to each file
        file1.write(','.join(column_titles) + '\n')
        file2.write(','.join(column_titles) + '\n')


        for core_sample, group in grouped:
            core_name = core_sample
            
            if core_name in groupings['Over or Equal to 60'] and (core_name in groupings['Grade 2'] or core_name in groupings['Grade 3']):
                group.to_csv('6K_RNA_Over60_Malignant.csv', mode = "a", index = False, header = False)

            if core_name in groupings['Under 60'] and (core_name in groupings['Grade 2'] or core_name in groupings['Grade 3']):
                group.to_csv('6K_RNA_Under60_Malignant.csv', mode = "a", index = False, header = False)

        # if core_name in groupings['Normal'] or core_name in groupings['NAT']:
        #     group.to_csv('6K_RNA_Normal_NAT.csv', mode = "a", index = False, header = False)

        # if core_name in groupings['Grade 2'] or core_name in groupings['Grade 3']:
        #     group.to_csv('6K_RNA_Malignant.csv', mode = "a", index = False, header = False)
        
        # if core_name in groupings['Basal']:
        #     group.to_csv('6K RNA Basal.csv', mode = "a", index = False, header = False)

        # if core_name in groupings['Classical']:
        #     group.to_csv('6K RNA Classical.csv', mode = "a", index = False, header = False)

        # if core_name in groupings['Over or Equal to 60']:
        #     group.to_csv('6K RNA Over or Equal to 60.csv', mode = "a", index = False, header = False)

        # if core_name in groupings['Under 60']:
        #     group.to_csv('6K RNA Under 60.csv', mode = "a", index = False, header = False)
        
