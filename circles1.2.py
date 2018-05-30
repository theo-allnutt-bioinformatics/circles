from __future__ import division
import sys
import os
import re
import glob
import copy
import math
from Tkinter import *
from tkFileDialog import *
from random import randint
from PIL import Image
from operator import itemgetter

#NB. blast outformat must be: -outfmt '6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore stitle'

#function to sort alphanumerically
digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

#get the bp range of a feature and return in a tuple				  
def getranges(l):
    r = []
    p = q = None
    for x in l + [-1]:
        if x - 1 == q:
            q += 1
        else:
            if p:
               if q > p:
                   r.append('%s-%s' % (p, q))
               else:
                   r.append('%s-%s' %(str(p),str(p)))
            p = q = x
    return '(%s)' % ', '.join(r)
	

#set the program options from the tkinter window				  
def get_options():

	if var  <>"":
		folder = var.get()
	else:
		print "no blast result folder specified, quitting"
		sys.exit(0)
		
	filelist=glob.glob(folder+"/*")
	filelist.sort(key=tokenize)
	print filelist
	
	if E2 <>"":
		winheight=int(E2.get())
	else:
		winheight = 5000
	if E1 <>"":
		winwidth=int(E1.get())
	else:
		winwidth=8000
	if E4 <>"":	
		refsize=int(E4.get())
	else:
		refsize=3036194
	if E5 <>"":	
		reftitle=E5.get()
	else:
		reftitle=""
	if E6 <>"":
		t2=float(E6.get())
	else:
		t2=70
	if E7 <>"":
		fontsize = int(E7.get())
	else:
		fontsize = 80
	if E8 <>"":	
		ann_file=E8.get()
	else:
		ann_file=""
	if E9 <>"":	
		thick=int(E9.get())
	else:
		thick=200	
	if E10 <>"":	
		dia=int(E10.get())
	else:
		dia=500
	if E11 <>"":	
		col=E11.get()
	else:
		col="235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30"
	
		
	draw_blast(filelist,winheight,winwidth,refsize,reftitle,t2,fontsize,offset,ann_file,thick,dia,anc,col)

#draw circles	
def draw_blast(filelist,winheight,winwidth,refsize,reftitle,t2,fontsize,offset,ann_file,thick,dia,anc,col):
	
	print "Drawing."
	print winwidth,winheight
	top2=Tk()
	bigwidth=0
	bigheight=0
	if winwidth>1400 or winheight>850: 
		bigwidth=winwidth
		bigheight=winheight
		winwidth=1400
		winheight=850
	
	if bigwidth<>0:
		frame=Frame(top2,width=winwidth,height=winheight)
		frame.grid(row=0,column=0)
		w1=Canvas(frame,height=int(winheight),width=int(winwidth),bg='white',scrollregion=(0,0,bigwidth,bigheight))
		hbar=Scrollbar(frame,orient=HORIZONTAL)
		hbar.pack(side=BOTTOM,fill=X)
		hbar.config(command=w1.xview)
		vbar=Scrollbar(frame,orient=VERTICAL)
		vbar.pack(side=RIGHT,fill=Y)
		vbar.config(command=w1.yview)
		w1.config(width=winwidth,height=winheight)
		w1.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
		w1.pack(side=LEFT,expand=True,fill=BOTH)
		winheight=bigheight
		winwidth=bigwidth
	else:	
		w1 = Canvas(top2,height=int(winheight),width=int(winwidth),confine=False,bg='white')
		w1.pack()
	
	colours=col.split(";")
	
	x0=(winwidth/2)-(winwidth*0.10)
	y0=winheight/2
	space=8
		####################################################### max ring thickness
	maxthick = ((winheight/2)-(dia+50)-(len(filelist)*space)-(fontsize*3))/len(filelist)
	print "max thickness=",maxthick
	if thick>maxthick:
		thick=maxthick
	
	dia0=dia+thick+space
	colours_used=[]
	c1=-1
	print filelist
	
	for f in filelist:
		span=0
		infile=open(f,'r')
		c1=c1+1
		
		#put file into dict
		data={}
		stpos=[]
		print 'reading and sorting file'
		for line in infile:
			pid = float(line.split("\t")[2])
			eval = float(line.split("\t")[9])
			start = int(line.split("\t")[7])
			end = int(line.split("\t")[8])
			span=span+(end-start)
			if start > end:	
				z= start
				start=end
				end=z
			span=span+(end-start)
			data[start,end]=[]
			data[start,end].append(start)
			data[start,end].append(end)
			data[start,end].append(pid)
			data[start,end].append(eval)
			stpos.append((start,end))
		print "Span",span	
		stpos=sorted(stpos,key=itemgetter(1),reverse=True)
		stpos=sorted(stpos,key=itemgetter(0))
		
		#delete internal hits
	
		#data2=copy.deepcopy(data)
		data2={}
		print "deleting internal hits"
		
		n = len(stpos)
		data2[stpos[0]]=data[stpos[0]]
		x=-1
		v2=1
		while x<n-2:
			x=x+v2
			v=x
			v2=1
			while x+v2<n:
				v=v+1
				if data[stpos[v]][1]<=data[stpos[x]][1]:
					v2=v2+1
				else:
					data2[stpos[v]]=data[stpos[v]] #b.append(a[v])
					break
			#print x,v,v2
		#get keys again
		stpos=[]
		
		for i in data2.keys():
			stpos.append(i)
			
		stpos=sorted(stpos,key=itemgetter(1),reverse=True)
		stpos=sorted(stpos,key=itemgetter(0)) #sort again
		
		#print stpos
		
		if c1>=len(colours)-1:
			print "Warning, not enough colours specified, using random."
			
			rc1=randint(20,235)
			rc2=randint(20,235)
			rc3=randint(20,235)
		
			randcol='%s,%s,%s' % (rc1, rc2, rc3)
			colours_used.append(randcol)
		else:
			colours_used.append(colours[c1])
		
		print "colour:",colours_used[c1]
		dia = dia+thick+space
		############scale parameter needed because pixels not square
		scale=12
		res1=int(scale*2*3.142*(dia+thick))  
		print f, res1
		gapdata=[]
		p1=-1
		print "Drawing hits"
		for k in stpos:
			p1=p1+1 #count position in key list stpos
			pid = float(data[k][2])
			eval = data[k][3]
			start = int(data[k][0])
			end = int(data[k][1])
			#populate gap file
			
			if p1>0:
				prevend=int(data[stpos[p1-1]][1])+1
				if start> prevend: #current start more than 1 more than previous end
					gapdata.append("%s-%s" %(str(prevend),str(start)))
			
			if pid>=float(t2):
				thick1=thick*pid/100
				#print pid,float(t2),start,end
				
				#need to make resolution the same on each circle
				#number of pixels in circle edge?
				st1=int(start/refsize*res1)
				en1=int(end/refsize*res1)
				for x in range(st1,en1+1):
					shad1=thick1*0.2
					shad2=thick1*0.8
					deg=(x/res1*360)-90
					rad=math.radians(deg)
					#
					x1=x0+(dia*math.cos(rad))
					y1=y0+(dia*math.sin(rad))
					x2=x0+((dia+thick1)*math.cos(rad))
					y2=y0+((dia+thick1)*math.sin(rad))
					
					xs2=x0+((dia+shad1)*math.cos(rad))
					ys2=y0+((dia+shad1)*math.sin(rad))
					
					xt2=x0+((dia+shad2)*math.cos(rad))
					yt2=y0+((dia+shad2)*math.sin(rad))
					
					if c1<=len(colours)-1:
						
						mc1=int(colours[c1].split(",")[0])
						mc2=int(colours[c1].split(",")[1])
						mc3=int(colours[c1].split(",")[2])
						mc='#%02x%02x%02x' % (mc1, mc2, mc3)
						
						shade=20
						ms1=mc1-shade
						ms2=mc2-shade
						ms3=mc3-shade
						ms='#%02x%02x%02x' % (ms1, ms2, ms3)
						
						mh1=mc1+shade
						mh2=mc2+shade
						mh3=mc3+shade
						mh='#%02x%02x%02x' % (mh1, mh2, mh3)
						
						rc1=int(colours_used[c1].split(",")[0])
						rc2=int(colours_used[c1].split(",")[1])
						rc3=int(colours_used[c1].split(",")[2])
						rc='#%02x%02x%02x' % (rc1, rc2, rc3)
						
						rs1=rc1-shade
						rs2=rc2-shade
						rs3=rc3-shade
						rs='#%02x%02x%02x' % (rs1, rs2, rs3)
						
						rh1=rc1+shade
						rh2=rc2+shade
						rh3=rc3+shade
						rh='#%02x%02x%02x' % (rh1, rh2, rh3)
						
						w1.create_line(x1,y1,x2,y2,fill=mc,width=0,state='disabled')
						#shadow
						w1.create_line(x1,y1,xs2,ys2,fill=ms,width=0,state='disabled')
						w1.create_line(xt2,yt2,x2,y2,fill=mh,width=0,state='disabled')
					else:
						w1.create_line(x1,y1,xs2,ys2,fill=rc,width=0,state='disabled')
						#shadow
						w1.create_line(x1,y1,xs2,ys2,fill=rs,width=0,state='disabled')
						w1.create_line(xt2,yt2,x2,y2,fill=rh,width=0,state='disabled')
					
		print "gaps=",len(gapdata)
		
		#Draw gaps
		print "Drawing gaps"
		
		bpt=refsize/((((thick+space)*len(filelist))+dia0)*3.142*2)
		print "gap threshold=", bpt, "bp"
		for i in gapdata:
			start = int(i.split("-")[0])
			end = int(i.split("-")[1])
			
			#print "gap resolution threshold=",bpt
			if end-start > bpt: #gapthreshold: #this is the threshold of pixels in bp
				st1=int(start/refsize*res1)
				en1=int(end/refsize*res1)
				for x in range(st1,en1+1):
					
					deg=(x/res1*360)-90
					rad=math.radians(deg)
						
					x1=x0+(dia*math.cos(rad))
					y1=y0+(dia*math.sin(rad))
					x2=x0+((dia+thick)*math.cos(rad))
					y2=y0+((dia+thick)*math.sin(rad))
					
					w1.create_line(x1,y1,x2,y2,fill='white',width=0,state='disabled')

		
	#add title
	w1.create_text(x0,y0-fontsize,text=reftitle,font=("Helvetica",fontsize))
	w1.create_text(x0,y0+(fontsize),text=str(refsize)+" bp",font=("Helvetica",fontsize-2))
	#add scale ring
	dia2=dia0-14
	w1.create_oval(x0-dia2,y0-dia2,x0+dia2,y0+dia2,width=3)
	d=-1
	dia3=dia2-15
	thick2=15
	for x in range(-90,269):
		deg=x
		rad=math.radians(deg)
		d=d+1
		bp1=int((deg+90)/360*refsize)
		
		
		if d==0 or d==45: #label every 45 degrees
			d=0
			x1=x0+(dia3*math.cos(rad))
			y1=y0+(dia3*math.sin(rad))
			x2=x0+((dia3+thick2)*math.cos(rad))
			y2=y0+((dia3+thick2)*math.sin(rad))	
			x3=x0+((dia3-18)*math.cos(rad))
			y3=y0+((dia3-18)*math.sin(rad))
			w1.create_line(x1,y1,x2,y2,fill="black",width=2)
			txt=str(int((deg+90)/360*refsize))
			if deg>225 or deg<-45:
				anc="n"
			if deg>=-45 and deg<=45:
				anc="e"
			if deg >45 and deg<135:
				anc="s"
			if deg >=135 and deg <=225:
				anc="w"
			w1.create_text(x3,y3,text=txt,font=("Helvetica",int(fontsize*0.65)),anchor=anc)
			
	#add annotations
	
	dia4=dia+thick+space #just outside outer ring
	dia5=dia4+(fontsize/2) #just before text, fontsize= LENGTH OF ANNOTATION LABEL LINE
	dia6=dia5+(fontsize/3) #text
	print "annotations",ann_file
	if ann_file<>"":
		print "annotating"
		print ann_file
		
		h=open(ann_file,'r')

		
		res1=int(scale*2*3.142*(dia4))
		
		for i in h:
			gene=i.split("\t")[0]
			coord=int(i.split("\t")[1].rstrip("\n"))
			print coord
			st1=int(coord/refsize*res1)
			deg1=(st1/res1*360)-90
			rad1=math.radians(deg1)
			
			x1=x0+(dia4*math.cos(rad1))
			y1=y0+(dia4*math.sin(rad1))
			x2=x0+(dia5*math.cos(rad1))
			y2=y0+(dia5*math.sin(rad1))
			x3=x0+(dia6*math.cos(rad1))
			y3=y0+(dia6*math.sin(rad1))
			
			if deg1>225 or deg1<=-45:
				anc="s"
			if deg1>-45 and deg1<=45:
				anc="w"
			if deg1 >45 and deg1<=135:
				anc="n"
			if deg1 >135 and deg1 <=225:
				anc="e"
			#Annotation text
			w1.create_text(x3,y3,text=gene,font=("Helvetica",int(fontsize*0.75)),anchor=anc)
			w1.create_line(x1,y1,x2,y2)
		dia6=dia6+30
	add_key(colours_used,dia6,filelist,winheight,winwidth,x0,y0,w1,fontsize)
	
	w1.bind("<Button-1>",lambda event, arg=[x0,y0,refsize,w1]: coords(event,arg))
	
def coords(event,arg):
	xh=int(arg[0])
	yh=int(arg[1])
	rf=arg[2]
	win=arg[3]
	
	x = win.canvasx(event.x)
	y = win.canvasy(event.y)
	
	if y==yh and x>xh: v=90
	if y>yh and x==xh: v=180
	if y==yh and x<xh: v=270
	if x==xh and y<yh: v=0
	
	if y <yh and x>xh: v = math.degrees(math.atan((x-xh)/(yh-y)))
	if y>yh and x>xh: v = 180-(math.degrees(math.atan((x-xh)/(y-yh))))
	if y>yh and x<xh: v = 180+(math.degrees(math.atan((xh-x)/(y-yh))))
	if y<yh and x<xh: v = 360-(math.degrees(math.atan((xh-x)/(yh-y))))
	
	print "%s bp" %(int(v/360*rf))
	
	if x>20 and y>20:
		win.create_polygon(18,8,18,32,88,32,88,8,fill = "white")
		win.create_text(20,20,anchor='w',font=("Helvetica",10),text="%s bp" %(int(v/360*rf)))
	else:
		win.create_polygon(18,8,18,32,88,32,88,8,fill = "white")
		
#add key
def add_key(colours_used,dia6,filelist,winheight,winwidth,x0,y0,w1,fontsize):
	print colours_used
	fontsize2 = fontsize*0.9
	defsize=fontsize2*1.3
	
	x=x0+dia6+(winwidth/25) #offset of key from right side of circle annotations
	
	n=len(colours_used)
	
	span = winheight*0.75
	
	divisions = span / n
	if divisions > defsize+(defsize*0.25):
		divisions = defsize+(defsize*0.25)
		
	gap = divisions*0.25
	
	boxy = divisions*0.75
	
	boxx = boxy*2
	
	keyfont = boxy/1.3
		
	start = y0-(span/2)

	y=start
	c=-1
	t5=0
	longest_name=0
	for i in colours_used:
		t5=t5+1
		c=c+1
		if len(filelist[c].split("\\")[-1].split(".")[0])>longest_name:
			longest_name=len(filelist[c].split("\\")[-1].split(".")[0])
		if t5>40:
			y=start
			x=x+(longest_name*1.3)+(boxx)+(keyfont*5)
			t5=0
		
		mk1=int(colours_used[c].split(",")[0])
		mk2=int(colours_used[c].split(",")[1])
		mk3=int(colours_used[c].split(",")[2])
		mk='#%02x%02x%02x' % (mk1, mk2, mk3)
		
		w1.create_polygon(x,y,x+boxx,y,x+boxx,y+boxy,x,y+boxy,fill = mk)
		w1.create_text(x+boxx+(keyfont+1.3),y+(boxy/2),anchor='w',font=("Helvetica",int(keyfont)),text=filelist[c].split("\\")[-1].split(".")[0])
		y=y+gap+boxy
	
	w1.update()
	
	i0=w1.winfo_rootx()
	j0=w1.winfo_rooty()
	i1=i0+w1.winfo_width()
	j1=j0+w1.winfo_height()
	
	w1.postscript(file="image.eps", colormode='color',width=winwidth,height=winheight)
	
	img= Image.open("image.eps")
	img.save("image.jpg",quality=100)
	#img.save("image.tiff")
	
# make a Tkinter menu for options
anc=""
top=Tk()

var = StringVar(top)
var.set("")
ann=StringVar(top)
ann.set("")
offset=100

L1 = Label(top, text="Window width")
L1.grid(row=1, column=1, sticky=W)
E1 = Entry(top)
E1.grid(row=1, column=2)

L2 = Label(top, text="Window height")
L2.grid(row=2, column=1, sticky=W)
E2 = Entry(top)
E2.grid(row=2, column=2)

L3 = Label(top, text="Blast input folder")
L3.grid(row=3, column=1, sticky=W)
E3 = Entry(top, textvariable=var)
E3.grid(row=3, column=2)
B1 = Button(top,text="Browse",command=lambda: var.set(askdirectory())).grid(row=3, column=3)

L4 = Label(top, text="Ref genome size")
L4.grid(row=4, column=1, sticky=W)
E4 = Entry(top)
E4.grid(row=4, column=2)

L5 = Label(top, text="Title")
L5.grid(row=5, column=1, sticky=W)
E5 = Entry(top)
E5.grid(row=5, column=2)

L6 = Label(top, text="%ID threshold")
L6.grid(row=6, column=1, sticky=W)
E6 = Entry(top)
E6.grid(row=6, column=2)

L7 = Label(top, text="Font size")
L7.grid(row=7, column=1, sticky=W)
E7 = Entry(top)
E7.grid(row=7, column=2)

L9 = Label(top, text="Ring thickness")
L9.grid(row=8, column=1, sticky=W)
E9 = Entry(top)
E9.grid(row=8, column=2)

L10 = Label(top, text="Inner radius")
L10.grid(row=9, column=1, sticky=W)
E10 = Entry(top)
E10.grid(row=9, column=2)

L8 = Label(top, text="Annotation file")
L8.grid(row=10, column=1, sticky=W)
E8 = Entry(top, textvariable=ann)
E8.grid(row=10, column=2)
B8=Button(top,text="Browse",command=lambda: ann.set(askopenfilename())).grid(row=10, column=3)

L11 = Label(top, text="Colour map")
L11.grid(row=11, column=1, sticky=W)
E11 = Entry(top)
E11.grid(row=11, column=2)


E1.insert(0,'8000')
E2.insert(0,'5000')
E3.insert(0,"C://python//test")
E4.insert(0,'3089739')
E5.insert(0,'ST204')
E6.insert(0,'50')
E7.insert(0,'80')
E8.insert(0,"C:\python\st204.ano")
E9.insert(0,'200')
E10.insert(0,'500')
E11.insert(0,"235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30")


B4 = Button(top,text="Draw", command=lambda: get_options()).grid(row=13, column=3)



raw_input()

#top.mainloop()


	

	

