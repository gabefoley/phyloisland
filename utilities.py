import os
import glob
def makeQueryString(iter, info = "", link = "", final = ""):
    """
    Take a list of items and concatanate them together
    :param iter: List of items to concatanate
    :param info: Additional piece of info to add behind each item
    :param link: Linking phrase to link all items
    :param final: Final phrase to appened to the enf of the string
    :return: String formatted with items and linking phrases
    """
    queryString = ""
    for item in iter:
        queryString += item + info + link

    # Remove the final joining string from the queryString
    queryString = queryString[:-len(link)] + final

    return queryString

def removeFile(*args):
    """
    Remove files in the list from the directory

    :param args: Files to remove
    :return:
    """
    for arg in args:
        os.remove(arg)