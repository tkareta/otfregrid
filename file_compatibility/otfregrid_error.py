class Error(Exception):
    # basing this off of Dreampy's
    # base error
    def __init__(self):
        pass

class otfregridfileerror(Error):

    def __init__(self, errorname, reason):
        # errorname: The name of the error for exception.
        # reason: An explanatory text that details the error
        # both of these are strings
 
        self.errorname = errorname
        self.reason = reason
        self.args = (self.reason,)
        self.message = errorname
        self.message += ": "
        self.message += reason
        self.str()
        
    def __str__(self):
        return self.message

class OTFRegridGeneralError(Error):

    def __init__(self, errorname, reason):
        # errorname: The name of the error for exception.
        # reason: An explanatory text that details the error
        # both of these are strings
 
        self.errorname = errorname
        self.reason = reason
        self.args = (self.reason,)
        self.message = errorname
        self.message += ": "
        self.message += reason
        self.str()
        
    def __str__(self):
        return self.message

class OTFRegridArgumentError(Error):

    def __init__(self, errorname, reason):
        # errorname: The name of the error for exception.
        # reason: An explanatory text that details the error
        # both of these are strings
 
        self.errorname = errorname
        self.reason = reason
        self.args = (self.reason,)
        self.message = errorname
        self.message += ": "
        self.message += reason
        self.str()
        
    def __str__(self):
        return self.message

