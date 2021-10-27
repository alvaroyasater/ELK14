from tkinter import *
from tkinter.filedialog import askopenfilename


def ViewFileName(filext="*"):
    """ View files - Choose one or cancel
        Input: filext = "*" - File extention in the first list
    """
    root = Tk()
    file = ""
    if filext == "*":
        fileclass = "All Files"
        fitype = "*.*"
    elif filext == "xlsx":
        fileclass = "Excel File"
        fitype = "*." + filext
    file = askopenfilename(filetypes=((fileclass, fitype),
                                      ("Text File", "*.txt*"),
                                      ("LP Files", "*.lp"),
                                      ("Excel Files", "*.xlsx")),
                           title="Choose a file")
    root.withdraw()
    print(file)
    return file
