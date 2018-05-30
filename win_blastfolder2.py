import sys
import re
import glob
import subprocess
from Tkinter import *
from tkFileDialog import *

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

# 'qseqid sseqid pident length qstart qend sstart send evalue bitscore'

def blastit(subj,fold,out,typ,top):
	print "Starting.."
	subj=subject.get()
	print subj
	fold=folder.get()
	#fold="C:/python/"
	print fold
	out=outputfolder.get()
	#out="C:/python/"
	print out
	typ=type.get()
	print typ
	
	#options=[subj,fold,out,typ,top]
	filelist=glob.glob(fold+"/*")

	filelist.sort(key=tokenize)
	print filelist
	for f in filelist:
		
		outname = out+"\\"+f.split("\\")[-1]+".blast"
		print outname
		if typ == "blastp":
			p1= subprocess.Popen('blastp -query %s -subject %s -out %s -outfmt "6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore" ' %(f,subj,outname), shell=True).wait() 

		if typ == "blastn":
			p1= subprocess.Popen('blastn -task blastn -query %s -subject %s -out %s -outfmt "6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore" ' %(f,subj,outname), shell=True).wait() 
		 
			
	print "Done.. blast again or ctrl-c to quit"	
		
top=Tk()
top.wm_title("Blast a folder utility")
subject = StringVar()
subject.set("reference")
folder=StringVar()
folder.set("C:/python/")
outputfolder = StringVar()
outputfolder.set("C:/python/")
type = StringVar()
type.set("blastn")

offset=100

L1 = Label(top, text="Reference fasta file")
L1.grid(row=1, column=1, sticky=W)
E1 = Entry(top, textvariable=subject)
E1.grid(row=1, column=2)
B1=Button(top,text="Browse",command=lambda: subject.set(askopenfilename())).grid(row=1, column=3)

L2 = Label(top, text="Folder of fasta files to blast")
L2.grid(row=2, column=1, sticky=W)
E2 = Entry(top, textvariable=folder)
E2.grid(row=2, column=2)
B2 = Button(top,text="Browse",command=lambda: folder.set(askdirectory())).grid(row=2, column=3)		

L3 = Label(top, text="Output folder")
L3.grid(row=3, column=1, sticky=W)
E3 = Entry(top, textvariable=outputfolder)
E3.grid(row=3, column=2)
B3 = Button(top,text="Browse",command=lambda: outputfolder.set(askdirectory())).grid(row=3, column=3)


L4 = Label(top, text="blast type: \nblastn, blastp, or blastx")
L4.grid(row=4, column=1, sticky=W)
E4 = Entry(top, textvariable=type)
E4.grid(row=4, column=2)

print "Starting.."
subj=subject.get()
print subj
fold=folder.get()
#fold="C://python//stx-phages"
print fold
out=outputfolder.get()
print out
typ=type.get()
print typ

options=[subj,fold,out,typ,top]

B4 = Button(top,text="Run Blast", command=lambda: blastit(*options)).grid(row=5, column=4)
'''
if subj and fold and out and typ:
	B4.config(state='normal')
else:
	B4.config(state='disabled')
'''
top.mainloop()




