def makeQueryString(iter, info = "", link = "", final = ""):
    queryString = ""
    for item in iter:
        queryString += item + info + link

    # Remove the final joining string from the queryString
    print ('This is the query string')
    print (queryString)

    queryString = queryString[:-len(link)] + final
    print ('This is the query string')
    print (queryString)
    return queryString