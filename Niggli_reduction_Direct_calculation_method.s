//------------------------------This algorithm is used to calculate Niggli reduced cell from a given primitive cell.------------------------------
//read the lattice parameters from the input string
void readcell(string initiocell,number &a0 ,number &b0 ,number &c0 ,number &d0 ,number &e0 ,number &f0)
{
number width=len(initiocell)	
number this,next,l,r,counter=0

string thischar,nextchar
for(number i=0;i<width-1;i++)
{
thischar=mid(initiocell,i,1)
nextchar=mid(initiocell,i+1,1)

if(thischar==" "||thischar==","||thischar==";"||thischar=="/"||thischar=="|")this=1
else this=0

if(nextchar==" "||nextchar==","||nextchar==";"||nextchar=="/"||nextchar=="|")next=1
else next=0

if(next-this==1)
{
Setpersistentnumbernote("ElectronDiffraction Tools:Temp:"+counter,i)
counter=counter+1
}

if(next-this==-1)
{
Setpersistentnumbernote("ElectronDiffraction Tools:Temp:"+counter,i)
counter=counter+1
}
}

string Firstchar=mid(initiocell,0,1)
if(Firstchar==" "||Firstchar==","||Firstchar==";"||Firstchar=="/"||Firstchar=="|")  
{
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:0",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:1",R)
a0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:2",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:3",R)
b0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:4",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:5",R)
c0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:6",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:7",R)
d0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:8",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:9",R)
e0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:10",L)
f0=val(right(initiocell,abs(width-L-1)))
}

else
{
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:0",R)
a0=val(left(initiocell,R+1))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:1",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:2",R)
b0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:3",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:4",R)
c0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:5",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:6",R)
d0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:7",L)
Getpersistentnumbernote("ElectronDiffraction Tools:Temp:8",R)
e0=val(mid(initiocell,l+1,abs(l-r)))

Getpersistentnumbernote("ElectronDiffraction Tools:Temp:9",L)
f0=val(right(initiocell,abs(width-L-1)))
}

deletepersistentnote("ElectronDiffraction Tools:Temp")
}


//Calculate angles between axsis
void  anglesBetweenAxis(number d0,number e0,number f0, number &nd, number &ne, number &nf)	//corresponding to α、β、γ
{
number d,e,f,n=0  //n is used to count the number of output.

//case 1: -1, 1, 1 for a, b, c    <==> 1,-1,-1  for a, b, c
d=d0
e=180-e0
f=180-f0

if((d<90&&e<90&&f<90)||(d>=90&&e>=90&&f>=90) ) //all acute or all obuse angle
{
nd=d
ne=e
nf=f

n=n+1
}

//case 2: 1, -1, 1 for a, b, c  <==> -1,1,-1  for a, b, c
d=180-d0
e=e0
f=180-f0

if((d<90&&e<90&&f<90)||(d>=90&&e>=90&&f>=90) ) //all acute or all obuse angle
{
nd=d
ne=e
nf=f

n=n+1
}

//case 3: 1, 1, -1 for a, b, c   <==> -1,-1,1  for a, b, c
d=180-d0
e=180-e0
f=f0

if((d<90&&e<90&&f<90)||(d>=90&&e>=90&&f>=90) ) //all acute or all obuse angle
{
nd=d
ne=e
nf=f

n=n+1
}

if(n==0) //other conditions
{
nd=d0
ne=e0
nf=f0
}
}
	
//------------------------------------The main program------------------------------------
//The input dialog
taggroup Inputdialog, Inputitems
Inputdialog=dlgcreatedialog("Direct calution of the Niggli Cell", Inputitems)
object dselsubdialog=alloc(UIFrame).init(Inputdialog)

string inputCell="9.0521, 6.3559, 12.0714, 163.00, 169.78, 19.82"
taggroup CellLabel=dlgcreatelabel("Input Cell")
taggroup Cell=dlgcreatestringfield(inputCell,50).dlgidentifier("InputCell")
taggroup CellGroup=dlggroupitems(CellLabel, Cell).dlgtablelayout(2,1,0)

taggroup NLabel=dlgcreatelabel("N")    //the range of indices u, v, w, 10 is enough to calculate Niggli cell
taggroup Nrange=dlgcreateIntegerfield(10,4).dlgidentifier("N")
taggroup NGroup=dlggroupitems(NLabel, Nrange).dlgtablelayout(2,1,0)

taggroup ParameterGroup=dlggroupitems(CellGroup, NGroup).dlgtablelayout(2,1,0)
Inputitems.dlgaddelement(ParameterGroup)

string initiocell
number epsval1,epsval2,epsval3, size
if(dselsubdialog.pose())
{
initiocell=Cell.DLGGetStringValue()
size=Nrange.dlggetvalue()	
}
else exit(0)

//Identify the input cell
number a0,b0,c0,d0,e0,f0  //the parameters of the input cell: a, b, c, α, β, γ
readcell(initiocell,a0,b0,c0,d0,e0,f0)

//define some parameters
number Vrec=a0*b0*c0*sqrt(1-(cos(Pi()/180*d0))**2-(cos(Pi()/180*e0))**2-(cos(Pi()/180*f0))**2+2*cos(Pi()/180*d0)*cos(Pi()/180*e0)*cos(Pi()/180*f0) )
number h=Vrec/(a0*b0*sin(Pi()/180*f0) )

result("\n---------------------------------------------------------- Niggli Reduction ----------------------------------------------------------\n")
result("The input cell: ("+a0+", "+b0+", "+c0+", "+d0+", "+e0+", "+f0+" | "+Vrec+"), Search Range: "+size+"\n")

//Define three vectors
number A1=a0  //Vector A: A=A1 x+A2 y+A3 z
number A2=0
number A3=0
number B1=b0*cos(Pi()/180*f0)   //Vector B
number B2=b0*sin(Pi()/180*f0)
number B3=0
number C1=c0*cos(Pi()/180*e0)   //Vector C
number C2=-c0*sin(Pi()/180*e0)*((cos(Pi()/180*e0)*cos(Pi()/180*f0)-cos(Pi()/180*d0) )/(sin(Pi()/180*e0)*sin(Pi()/180*f0) ) )   //用模计算y方向投影
number C3=h

number a,b,c,d,e,f
a=a0**2  
b=b0**2  
c=c0**2  

d=2*b0*c0*cos(Pi()/180*d0) 
e=2*a0*c0*cos(Pi()/180*e0)  
f=2*a0*b0*cos(Pi()/180*f0)  

Deletepersistentnote("ElectronDiffraction Tools:Temp")  //Clear the tags

number t,u,v,w,i,j,xsize,ysize
xsize=(2*size+1) //u, v, w include negative, zero and position.
ysize=2*size +1  //u, v, w include negative, zero and position.

image imguv=realimage("",4,xsize*xsize,ysize)  //save data u, v, t
image imgW=imguv*0   ////save data u, v, w
number z0=trunc(xsize*ysize/2)

//-------------Calculate the size of the translation vectors-------------
for(w=-size;w<=size;w++)
{
for(v=-size;v<=size;v++)
{
for(u=-size;u<=size;u++)
{
t=sqrt(u**2*a+v**2*b+w**2*c+v*w*d+u*w*e+u*v*f)

i=z0+w*xsize+u
j=size+v

if(u==0&&v==0&&w==0)
{
imguv.SetPixel(i,j,1e6)  //Keep 0 0 0 is maximum
imgW.SetPixel(i,j,w)
}

else
{
imguv.SetPixel(i,j,t)  //Save u, v, w, and t.
imgW.SetPixel(i,j,w)
}
}
}
}

//-------------------------And now we will search three non-coplanar minimum edges---------------------------
image IndiciesMatrix=realimage("",4,3,3)  //save hkl indices
number u1,v1,w1,u2,v2,w2,u3,v3,w3,t1,t2,t3,x,y,dist
number A1r,A2r,A3r,B1r,B2r,B3r,C1r,C2r,C3r  // coordinates of the new vectors
number AngleD,AngleE,AngleF

//Define three images to save vector parameters
image AVector=realimage("",4,3,1)
image BVector=realimage("",4,3,1)
image CVector=realimage("",4,3,1)

number no1,no2,no3
number n1=1,n2=1,n3=1
while(n1>0) 
{
//find the first minimum value of t: t1
t1=imguv.min(x,y)
imguv.SetPixel(x,y, imguv.max())   

w1=imgW.GetPixel(x,y)
u1=x-z0-w1*xsize
v1=y-size
imguv.SetPixel(z0-w1*xsize-u1,size-v1, imguv.max())   // Note -u, -v, -w has the same t with u, v, w, skip it

while(n2>0)
{
//find the second minimum value of t: t2
t2=imguv.min(x,y)
imguv.SetPixel(x,y, imguv.max())   

w2=imgW.GetPixel(x,y)
u2=x-z0-w2*xsize
v2=y-size
imguv.SetPixel(z0-w2*xsize-u2,size-v2, imguv.max())   // Note -u, -v, -w has the same t with u, v, w, skip it

IndiciesMatrix.SetPixel(0,0,u1)  //save hkl indices
IndiciesMatrix.SetPixel(0,1,v1)
IndiciesMatrix.SetPixel(0,2,w1)

IndiciesMatrix.SetPixel(1,0,u2)
IndiciesMatrix.SetPixel(1,1,v2)
IndiciesMatrix.SetPixel(1,2,w2)

//The new vectors OA and OB
A1r=u1*A1+v1*B1+w1*C1  
A2r=u1*A2+v1*B2+w1*C2
A3r=u1*A3+v1*B3+w1*C3

B1r=u2*A1+v2*B1+w2*C1  
B2r=u2*A2+v2*B2+w2*C2
B3r=u2*A3+v2*B3+w2*C3

AVector.setpixel(0,0,A1r)
AVector.setpixel(1,0,A2r)
AVector.setpixel(2,0,A3r)

BVector.setpixel(0,0,B1r)
BVector.setpixel(1,0,B2r)
BVector.setpixel(2,0,B3r)

image AB=CrossProduct(AVector,BVector)
if(AB.GetPixel(0,0)!=0||AB.GetPixel(1,0)!=0||AB.GetPixel(2,0)!=0)   //Make sure  two vectors t1 and t2 are non-collinear，010,020,0-20; 111,222
{
while(n3>0)
{
//find the third minimum value of t: t3
t3=imguv.min(x,y)
imguv.SetPixel(x,y, imguv.max())   

w3=imgW.GetPixel(x,y)
u3=x-z0-w3*xsize
v3=y-size
imguv.SetPixel(z0-w3*xsize-u3,size-v3, imguv.max())   //Note -u, -v, -w has the same t with u, v, w, skip it

IndiciesMatrix.SetPixel(2,0,u3)
IndiciesMatrix.SetPixel(2,1,v3)
IndiciesMatrix.SetPixel(2,2,w3)

C1r=u3*A1+v3*B1+w3*C1  //the third vector
C2r=u3*A2+v3*B2+w3*C2
C3r=u3*A3+v3*B3+w3*C3

CVector.setpixel(0,0,C1r)
CVector.setpixel(1,0,C2r)
CVector.setpixel(2,0,C3r)

number Na,Nb,Nc,Nd,Ne,Nf,V
Na=t1
Nb=t2
Nc=t3

AngleD=acos(DotProduct(BVector,CVector)/(t2*t3))*180/Pi()
AngleE=acos(DotProduct(AVector,CVector)/(t1*t3))*180/Pi()
AngleF=acos(DotProduct(AVector,BVector)/(t1*t2))*180/Pi()

anglesBetweenAxis(AngleD,AngleE,AngleF, Nd,Ne,Nf)
V=abs(DotProduct(AVector,CrossProduct(BVector,CVector)))
//V=na*nb*nc*sqrt(1-(cos(nd*Pi()/180))**2-(cos(ne*Pi()/180))**2-(cos(nf*Pi()/180))**2+2*cos(nd*Pi()/180)*cos(ne*Pi()/180)*cos(nf*Pi()/180) )

if(V>1e-6)
{
result("1st: "+format(t1, "%8.4f" )+" ("+format(u1, "%3.0f" )+", "+format(v1, "%3.0f" )+", "+format(w1, "%3.0f" )+"); ("+format(A1r, "%8.4f" )+", "+format(A2r, "%8.4f" )+", "+format(A3r, "%8.4f" )+")\n")
result("2nd: "+format(t2, "%8.4f" )+" ("+format(u2, "%3.0f" )+", "+format(v2, "%3.0f" )+", "+format(w2, "%3.0f" )+"); ("+format(B1r, "%8.4f" )+", "+format(B2r, "%8.4f" )+", "+format(B3r, "%8.4f" )+")\n")
result("3th: "+format(t3, "%8.4f" )+" ("+format(u3, "%3.0f" )+", "+format(v3, "%3.0f" )+", "+format(w3, "%3.0f" )+"); ("+format(C1r, "%8.4f" )+", "+format(C2r, "%8.4f" )+", "+format(C3r, "%8.4f" )+")\t")
result(" ==> "+format(Na, "%8.4f" )+", "+format(Nb, "%8.4f" )+", "+format(Nc, "%8.4f" )+" | "+format(Nd, "%8.2f" )+", "+format(Ne, "%8.2f" )+", "+format(Nf, "%8.2f" )+" | "+format(V, "%8.4f" )+"\n\n")

n2=0
n1=0
break
}
}
}
}
}
