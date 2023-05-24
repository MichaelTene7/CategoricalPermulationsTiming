## USAGE:
# accepts plain text or code command line arguments, and returns the argument's result as a string. 
#The argument should be in the format of: x=ARGVAL
# "marker" is the letter used to indicate which argument is which. in the above example, the marker is <x>, in the below example, it is <j>
# "message" is the message sent if the argument is missing


# example:
# Command line argument: j=MyString   Or   j=<CodeToEvaluateHere>
# jValue = cmdArgImport('j', message = "my message for a missing argument here")
# Returns: jValue = MyString        Or       jValue = <EvaluatedCode>

cmdArgImport = function(marker, message = ""){
  marker                                                           #send the marker value
  markerWhole = paste("^", marker, "=", sep='')                    #convert marker to grep format
  commandLineValue = grep(markerWhole, args, value = TRUE)         #get a string based on the identifier   
  if(length(commandLineValue != 0)){                               #If the string is not empty:
    commandLineString = substring(commandLineValue, 3)             #make a string without the identifier
    if(grepl('(', commandLineValue, fixed = TRUE)){                # if the input is code -- has a '(' in it
      commandLineOutput = eval(str2lang(commandLineString))        #convert that string to code, then evaluate that code
    }else{                                                         #else
      commandLineOutput = commandLineString                        #use the string directly 
    }
    message(commandLineOutput)                                     #Report the result
    return(commandLineOutput)
  }else{                                                           #if it is empty, send a message; return NA
    if(message == ""){                                             #if mo message specified
      defaultMessage = paste("Arg associated with [", marker, "] not specified") # set default message
      message = defaultMessage                                     #use default message
      }                    
    message(message)
    return(NA)
  }
}

#testOut = cmdArgImport('n')
