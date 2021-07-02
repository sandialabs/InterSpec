from ctypes import *
import os
import time
import urllib
import random
#import urllib.parse #python3
import webbrowser

# This file demonstrates starting an InterSpec server instance from Python, and two ways to load
#  spectrum files from the filesystem.
# Not tested beyond "seems to niavely, kinda, sorta work", good luck


cwd = os.getcwd()
sharedLibFileName = os.path.join(cwd,"InterSpec.dylib")
exampleFileName = os.path.join(cwd,"example.n42");
interspec = CDLL( sharedLibFileName )


process_name = "InterSpec"
user_data_dir = "userdata"
basedir = "."
xml_config_path = "data/config/wt_config_electron.xml"
port = interspec.interspec_start_server( process_name, user_data_dir, basedir, xml_config_path )

# First method to open a file: just include the file path in the URL
url = "http://localhost:" + str(port) + "?" + urllib.urlencode({'specfilename': exampleFileName})
#url = "http://localhost:" + str(port) + "?" + urllib.parse.urlencode({'specfilename': exampleFileName})  #python3
print( "pointing your browser to " + url )

# On my computer at least, there seems to eb a built-in delay for opening new tabs in the browser
webbrowser.open(url)

time.sleep(15)

print( "Will now open a second session that will be manipulated after opening" ) 

# Second method to open a file: create a session with a unique token, and then use that token to to manipulate 
#  that session from the code here.  Doesnt matter what this token string is, as long as it is unique
#  per running session.
token = str(random.getrandbits(64))
interspec.interspec_add_allowed_session_token(token)

# Because we are re-using the code written to interface to Electron , the "primary=no" argument
#  is needed to keep the app from trying to create native Electron menus (instead html-based
#  menus will be created)
url = "http://localhost:" + str(port) + "?restore=no&primary=no&apptoken=" + token
webbrowser.open(url)

# You will see a bunch of debug messages printout from the c++ on the terminal here - sorry!

# Wait for the new session to load and stuff
while interspec.interspec_session_is_alive(token) == 0:
    time.sleep(1)
print("Your session has been loaded in the browser")

print( "Will wait 10 seconds for you to realize the new empty session" )
time.sleep(10)


# We can open multiple files up at a time (although this doesnt really make sense here),
#  thus this next call expects the input to be a JSON array of strings that give paths
#  to the file to open 
fileOpenJson = '["' + exampleFileName + '"]'

print( "Will open example.n42" )
interspec.interspec_open_file( token, fileOpenJson )

# You'll see a bunch more debug printout from C++ again here

# Right now loading files is the only action defined for interacting with specific 
#  sessions from the code here, but it would be pretty easy to add more actions.

print( "You have one minute - better analyze quickly!" ) 
time.sleep(60)

print( "Killing server" )
interspec.interspec_kill_server()

print( "exiting" )
