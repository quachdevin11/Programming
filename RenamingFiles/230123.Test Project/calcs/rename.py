import os

# Gets folder path and name of the project
folder_path = os.getcwd()
name_cwd = folder_path.split("\\")[-2]

# Type of files in Calcs folder to rename 
filetypes_to_rename = ['Cover', 'Letter', 'GlazingCalcs']

# Steps through each file in calcs folder, renames the job number and project name of the type of files above. Ex. 123456.Test Project.Cover.docx --> type = Cover (Always second from last)
for file in os.listdir(folder_path):
    number = file.split('.')[0]
    type = file.split('.')[1]
    if type in filetypes_to_rename:
        old_name = os.path.join(folder_path, file)
        new_name = old_name.replace(number, name_cwd)
        os.rename(old_name, new_name)

