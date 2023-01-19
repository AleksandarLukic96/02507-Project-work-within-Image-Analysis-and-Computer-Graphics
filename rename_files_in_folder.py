# importing os module
import os

path = os.getcwd()
 
# Function to rename multiple files
def main():   
    folder = path + "\\gif\\"
    for count, filename in enumerate(os.listdir(folder)):
        dst = f"SP_Frame_{str(count)}.png"
        src =f"{folder}/{filename}"  # foldername/filename, if .py file is outside folder
        dst =f"{folder}/{dst}"
         
        # rename() function will rename all the files
        os.rename(src, dst)
 
# Driver Code
if __name__ == '__main__':
     
    # Calling main() function
    main()