
import sys
import os
import subprocess

def bash( bashCommand ):
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        #process = subprocess.Popen(bashCommand.split())
        output, error = process.communicate()
        return output ,error

path = sys.argv[1]
classname = sys.argv[2]



out = bash("ls "+path)
out = out[0].split("\n")
out = [ x for x in out if '.root' in x]


f = open(classname+".C","w")

f.write("class "+classname+"{\n")
f.write("	public:\n")
f.write("	std::string PATH=\""+path+"\";\n")
f.write("	std::vector<std::string> files{\n")
f.write("		")
for i in range(len(out)-1):
	f.write("\""+out[i]+"\", ")
f.write("\""+out[-1]+"\" };\n")
f.write("};\n")
