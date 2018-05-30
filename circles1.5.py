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
	
	
	#sort filelist by total blast score...
	filebits=[]
	for i in filelist:
		bits=0
		f5=open(i,'r')
		for x in f5:
			#print x
			bits=bits+float(x.split("\t")[-2])
		f5.close()	
		filebits.append((i,bits))
	
	bitsorted=sorted(filebits,key=itemgetter(1),reverse=True)
	filelist2=[]
	for i in bitsorted:
		filelist2.append(i[0])
	
	
	
	if E2 <>"":
		winheight=int(E2.get())
	else:
		winheight = 900
	if E1 <>"":
		winwidth=int(E1.get())
	else:
		winwidth=1100
	if E4 <>"":	
		ref_end=int(E4.get())
	else:
		ref_end=1000000
	if E5 <>"":	
		reftitle=E5.get()
	else:
		reftitle=""
	if E6 <>"":
		t2=float(E6.get())
	else:
		t2=50
	if E7 <>"":
		ref_ann = E7.get()
	else:
		ref_ann=""
	if E8 <>"":	
		ann_file=E8.get()
	else:
		ann_file=""
	if E9 <>"":	
		ref_start=int(E9.get())
	else:
		ref_start=0	
	if E10 <>"":	
		c_file=E10.get()
		col_file=open(c_file,'r')
		col=col_file.readline()
	else:
		col="235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30,235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30,235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30"
	
		
	draw_blast(filelist,filelist2,winheight,winwidth,ref_start,ref_end,reftitle,t2,offset,ann_file,ref_ann,anc,col)

#draw circles	
def draw_blast(filelist,filelist2,winheight,winwidth,ref_start,ref_end,reftitle,t2,offset,ann_file,ref_ann,anc,col):
	
	refsize=ref_end-ref_start
	
	#get annotation data
	if ref_ann<>"":
		f1=open(ref_ann,'r')
		
		global ann_data
		
		
		ann_data={}
		for i in f1:
			
			if i.split("\t")[0]=="CDS":
				cds_start = int(i.split("\t")[7])
				cds_end = int(i.split("\t")[8])
				gene = i.split("\t")[13]#+":"+i.split("\t")[15]
				
				for x in range(cds_start,cds_end):
					ann_data[x]=gene
														
	
	print "Drawing."

	top2=Tk()
	top2.title("Circles alignment viewer")
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
	
	space=int(winheight/900) # gap between rings
	if space <1:
		space=1
	print 'space',space
	#if fontsize>winheight/60:
	fontsize=int(winheight/55)
	#if dia>winheight/6: #set max inside ring radius
	
	dia=int(winheight/7)
	print 'radius', dia
	
	x0=(winwidth/2.25)#left offset so that key can fit on the right
	y0=winheight/2	
		####################################################### max ring thickness
	n_rings = len(filelist)
	if n_rings <5:
		n_rings =5
	maxthick = ((winheight/2)-(dia+20)-(len(filelist)*space)-(fontsize*3))/n_rings
	print "max thickness=",maxthick
	#if thick>maxthick:
	thick=maxthick
	
	dia0=dia #innermost radius
	
	colours_used=[]
	c1=-1
	#print filelist2
	
	#draw rings
	dia = dia+(thick/2)+space
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
			#span=span+(end-start)
			if start > end:	
				z= start
				start=end
				end=z
			span=span+(end-start)
		################## collect data from specified region- including hits that overlap it	
			
			if start <= ref_start and end >= ref_start and end <=ref_end:
				start=ref_start
				#end=end
			
			if start <= ref_start and end >= ref_end:
				end = ref_end
				start = ref_start
			
			if start >= ref_start and end >= ref_end:
				end = ref_end
				#start=start
			
			if start >= ref_start and start <=ref_end and end >= ref_start:
			
				data[start,end]=[]
				data[start,end].append(start)
				data[start,end].append(end)
				data[start,end].append(pid)
				data[start,end].append(eval)
				stpos.append((start,end))	
				
		#sort by reverse end pos then by start pos	
		stpos=sorted(stpos,key=itemgetter(1),reverse=True)
		stpos=sorted(stpos,key=itemgetter(0))
		#print stpos
		
		#delete internal hits
	
		#data2=copy.deepcopy(data)
		data2={}
		#print "deleting internal hits - blast hits that are within larger ones and therefore should be ignored"
		#print stpos 
		
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
		stpos2=[]
		for i in data2.keys():
			stpos.append(i)
			
		stpos=sorted(stpos,key=itemgetter(1),reverse=True)
		stpos=sorted(stpos,key=itemgetter(0)) #sort again
		
		#remove position data overlaps
		data3={}
		stpos2.append(stpos[0])
		for i in range(1,len(stpos)):
			if stpos[i][0]<stpos[i-1][1]:
				
				stpos2.append((stpos[i-1][1],stpos[i][1]))
			else:
				stpos2.append(stpos[i])
			
		c9=-1
		for i in stpos2:
			c9=c9+1
			data3[i]=data2[stpos[c9]]
			
		#print 'internal hits removed',stpos2
		#print data3
		
		
		
		###############get ring colours
		if c1>=len(colours)-1:
			print "Warning, not enough colours specified, using random."
			
			rc1=randint(20,235)
			rc2=randint(20,235)
			rc3=randint(20,235)
		
			randcol='%s,%s,%s' % (rc1, rc2, rc3)
			colours_used.append(randcol)
		else:
			colours_used.append(colours[c1])
		
		#print "colour:",colours_used[c1]
		
		if c1<=len(colours)-1:
			
			shade=20
			
			mc1=int(colours[c1].split(",")[0])
			mc2=int(colours[c1].split(",")[1])
			mc3=int(colours[c1].split(",")[2])
			
			if mc1-shade<=shade:mc1=mc1+shade
			if mc2-shade<=shade:mc2=mc2+shade
			if mc3-shade<=shade:mc3=mc3+shade
			if mc1+shade>=255:mc1=mc1-shade
			if mc2+shade>=255:mc2=mc2-shade
			if mc3+shade>=255:mc3=mc3-shade
			
			mc='#%02x%02x%02x' % (mc1, mc2, mc3)
			
			ms1=mc1-shade
			ms2=mc2-shade
			ms3=mc3-shade
			ms='#%02x%02x%02x' % (ms1, ms2, ms3)
				
			mh1=mc1+shade
			mh2=mc2+shade
			mh3=mc3+shade
			mh='#%02x%02x%02x' % (mh1, mh2, mh3)
			
			#print 'rgb',mc1,mc2,mc3
			#print 'hex',mc
			
		else: #not enough colours specified - choose random ones for rings 
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
		
			mc = rc
			ms = rs
			mh=rh

		scale=1
		res1=2*3.142*dia #circumference in pixels at middle
		resbp=refsize/res1
		print f
		print 'resolution=',resbp,'bp'
		print 'radius pixels',dia
		print 'circumference pixels', res1
		print 'thickness pixels', thick
		gapdata=[]
		p1=-1
		#print "Drawing hits"
		data={}
		data2={}
		for k in stpos2:
			p1=p1+1 #count position in key list stpos
			pid = float(data3[k][2])
			eval = data3[k][3]
			start = k[0]
			end = k[1]
			#populate gap file
			
			if p1>0:
				prevend=stpos2[p1-1][1]+1
				if start> prevend: #current start more than 1 more than previous end
					gapdata.append("%d-%d" %(prevend,start))
			'''
			#make gaps large enough to show#####################
			if p1>0:
				prevend=stpos2[p1-1][1]+1
				if start> prevend: #current start more than 1 more than previous end	
					start=start + int(resbp)
			'''
			st1=(start/refsize)
			en1=(end/refsize)	
					
			#print 'pid',pid,'t2',float(t2)
			if pid>=float(t2):
			
				thick1=thick*(pid/100) #n.b. arc drawing thickness is from centre of arc
				res1=2*3.142*dia#circumference in pixels middle
				res2=2*3.142*(dia-(thick1/2))#circumference in pixels inner edge
				res3=2*3.142*(dia+(thick1/2))#circumference in pixels outer edge
				'''
				#remove glitches in tkinter arc that makes every 90 degree arc have rounded ends
				#print st1,st1*360, refsize/res1
				if st1*360>(0-(refsize/res1)) and st1*360<(0+(refsize/res1)):
					st1=(st1+(refsize/res1))/refsize
					#print 'fixed',st1
				'''
				arcstart=(360-(st1*360))+90
				arcend=(360-(en1*360))+90
				
				arcextent=arcend-arcstart
			
				#deg=(x/res1*360)-90
				#rad=math.radians(deg)
			
			#calc arc parameters
				xarc0=x0-dia#-(thick1/2) #top left
				yarc0=y0-dia#-(thick1/2)
				xarc1=x0+dia#+(thick1/2) #bottom right
				yarc1=y0+dia#+(thick1/2)
				
				shadow=thick1*0.1
				#inner shadow
				xs0=xarc0+(thick1/2)-(shadow/2)
				ys0=yarc0+(thick1/2)-(shadow/2)
				xs1=xarc1-(thick1/2)+(shadow/2)
				ys1=yarc1-(thick1/2)+(shadow/2)
				
				#outer shadow
				xt0=xarc0-(thick1/2)+(shadow/2)
				yt0=yarc0-(thick1/2)+(shadow/2)
				xt1=xarc1+(thick1/2)-(shadow/2)
				yt1=yarc1+(thick1/2)-(shadow/2)
				
				#print arcstart, arcend, -arcextent
				
				#print ms,mc,mh
				
				if (-arcextent/360)*res2>1: #test if arc is longer than one pixel
					w1.create_arc(xarc0,yarc0,xarc1,yarc1,start=arcstart,extent=arcextent,width=thick1,outline=mc,style="arc",state='disabled')
				#inner shadow
					w1.create_arc(xs0,ys0,xs1,ys1,start=arcstart,extent=arcextent,width=shadow,outline=ms,style="arc",state='disabled')
				#outer shadow
					w1.create_arc(xt0,yt0,xt1,yt1,start=arcstart,extent=arcextent,width=shadow,outline=mh,style="arc",state='disabled')
				#w1.update()
				
		#print "gaps=",len(gapdata)
		

				#Draw gaps
		
		print "Drawing gaps"
		
		bpt=resbp#refsize/((((thick+space)*len(filelist))+dia0)*3.142*2)
		print "gap threshold=", bpt, "bp"
		print len(gapdata),"gaps"
		for i in gapdata:
			startg = int(i.split("-")[0])
			endg = int(i.split("-")[1])
			
			if startg >=ref_start and endg<= ref_end:
				
				#print "gap resolution threshold=",bpt
				if endg-startg > bpt: #gapthreshold: #this is the threshold of pixels in bp
					st1=int(startg/refsize*res1)
					en1=int(endg/refsize*res1)
					for x in range(st1,en1+1):
						
						deg=(x/res1*360)-90
						rad=math.radians(deg)
							
						x1=x0+((dia-(thick/2))*math.cos(rad))
						y1=y0+((dia-(thick/2))*math.sin(rad))
						x2=x0+((dia+(thick/2))*math.cos(rad))
						y2=y0+((dia+(thick/2))*math.sin(rad))
						
						w1.create_line(x1,y1,x2,y2,fill='white',width=0)
		
		dia = dia+thick+space ####################increment diameter for each ring
		
	#add title
	w1.create_text(x0,y0-int(fontsize),text=reftitle,font=("Helvetica",int(fontsize*0.8)))
	w1.create_text(x0,y0+int(fontsize),text=str(str(ref_start)+"-"+str(ref_end))+" bp",font=("Helvetica",int(fontsize*0.8)))
	
	#add scale ring
	dia2=dia0-(space*2)
	
	w1.create_oval(x0-dia2,y0-dia2,x0+dia2,y0+dia2,width=3)
	
	d=-1
	dia3=dia2-(space*3)
	
	for x in range(-90,269):
		deg=x
		rad=math.radians(deg)
		d=d+1
		bp1=int((deg+90)/360*refsize)
		
		if d==0 or d==45: #label every 45 degrees
			d=0
			#dash mark
			x1=x0+(dia3*math.cos(rad))
			y1=y0+(dia3*math.sin(rad))
			x2=x0+((dia2)*math.cos(rad))
			y2=y0+((dia2)*math.sin(rad))	
			
			#bp scale ring
			x3=x0+((dia3-(space*2))*math.cos(rad))
			y3=y0+((dia3-(space*2))*math.sin(rad))
			
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
			w1.create_text(x3,y3,text=txt,font=("Helvetica",int(fontsize*0.55)),anchor=anc)
			
	#add annotations to outer ring
	
	dia4=dia #just outside outer ring
	dia5=dia4+thick #just before text, fontsize= LENGTH OF ANNOTATION LABEL LINE
	dia6=dia5+space #text
	#print "annotations",ann_file
	if ann_file<>"":
		#print "annotating"
		#print ann_file
		
		h=open(ann_file,'r')

		res1=int(scale*3.142*(dia4*2))
		
		for i in h:
			gene=i.split("\t")[0]
			coord=int(i.split("\t")[1].rstrip("\n"))
			if coord>=ref_start and coord<=ref_end:
				if len(i.split("\t"))>2:
					ancol=i.split("\t")[2]
				else:
					ancol="Black"
				#print coord
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
				w1.create_text(x3,y3,text=gene,font=("Helvetica",int(fontsize*0.8)),anchor=anc,fill=ancol)
				w1.create_line(x1,y1,x2,y2,width=2)
		dia6=dia6+30
	add_key(colours_used,dia6,filelist,filelist2,winheight,winwidth,x0,y0,w1,fontsize)
	
	top3=Tk()
	top3.title("Annotations")
	popup = Canvas(top3,height=100,width=250,confine=False,bg='white')
	popup.pack()
	
	w1.bind("<Button-1>",lambda event, arg=[x0,y0,refsize,w1,ref_start]: get_coords(event,arg,fontsize,refsize,winwidth,popup))

def get_coords(event,arg,fontsize,refsize,winwidth,popup):
	xh=int(arg[0])
	yh=int(arg[1])
	rf=arg[2]
	win=arg[3]
	ref_start=arg[4]
	
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
	
	coords=int(v/360*rf)+ref_start
	
	print "%s bp" %(str(coords)), ann_data[coords]
	
	annfont=12
	
	popup.delete("all")
	
	if x>20 and y>20:

		popup.create_text(18,8,anchor='nw',font=("Helvetica",annfont),text=str(coords)+" bp: "+ann_data[coords],width=220)
		
		
	
		
#add key
def add_key(colours_used,dia6,filelist,filelist2,winheight,winwidth,x0,y0,w1,fontsize):
	#print colours_used
	fontsize2 = fontsize*0.9
	defsize=fontsize2*1.3
	
	x=x0+dia6+(winwidth/18) #offset of key from right side of circle annotations
	
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

top.title("Circles v1.4")
#define variables
var = StringVar(top)
var.set("")
ann=StringVar(top)
ann.set("")
col=StringVar(top)
col.set("")
ref_ann=StringVar(top)
ref_ann.set("")

offset=100
#global ann_data

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

L4 = Label(top, text="Ref genome End")
L4.grid(row=5, column=1, sticky=W)
E4 = Entry(top)
E4.grid(row=5, column=2)

L5 = Label(top, text="Title")
L5.grid(row=6, column=1, sticky=W)
E5 = Entry(top)
E5.grid(row=6, column=2)

L6 = Label(top, text="%ID threshold")
L6.grid(row=7, column=1, sticky=W)
E6 = Entry(top)
E6.grid(row=7, column=2)

L7 = Label(top, text="Ref annotation file")
L7.grid(row=8, column=1, sticky=W)
E7 = Entry(top, textvariable=ref_ann)
E7.grid(row=8, column=2)
B7=Button(top,text="Browse",command=lambda: ref_ann.set(askopenfilename())).grid(row=8, column=3)

L8 = Label(top, text="Outer ring annotation file")
L8.grid(row=9, column=1, sticky=W)
E8 = Entry(top, textvariable=ann)
E8.grid(row=9, column=2)
B8=Button(top,text="Browse",command=lambda: ann.set(askopenfilename())).grid(row=9, column=3)

L9 = Label(top, text="Ref genome start")
L9.grid(row=4, column=1, sticky=W)
E9 = Entry(top)
E9.grid(row=4, column=2)

L10 = Label(top, text="Colour map file")
L10.grid(row=10, column=1, sticky=W)
E10 = Entry(top, textvariable=col)
E10.grid(row=10, column=2)
B10=Button(top,text="Browse",command=lambda: ann.set(askopenfilename())).grid(row=10, column=3)

E1.insert(0,'8500')
E2.insert(0,'5000')
E3.insert(0,"C:\\python\\vibrio_blast")
E4.insert(0,'5152433')
E5.insert(0,'V. parahaemolyticus ATCC_17802')
E6.insert(0,'50')
E7.insert(0,"")
E8.insert(0,"")
E9.insert(0,'0')
E10.insert(0,"C:\\python\\ming_colormap.txt")

#colour_map="235,128,114;46,139,87;235,193,203;235,215,20;20,20,128;221,160,221;106,90,205;95,158,160;234,164,96;165,42,42;30,144,232;152,231,152;235,130,235;228,127,80;128,128,20;235,228,225;173,235,47;211,211,211;205,92,92;20,235,127;235,235,20;235,20,235;210,105,30"

B4 = Button(top,text="Draw", command=lambda: get_options()).grid(row=11, column=2)



#raw_input()

top.mainloop()


	

	

