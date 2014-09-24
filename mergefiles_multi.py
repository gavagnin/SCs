#######################################################################################
#
#                                  MERGEFILES_MULTI 
#                   Initial Condition creator for SCs merger simulation
#
# files to be merged= ftbm
# calling sintax: mergefiles_multi.py #ftbm name_ftbm name_output ICconfiguration extent
#
#ICconfiguration: sphere or filament
#extent= radius of the sphere or length of filament 
#
########################################################################################
import os
import gzip
import string
import re   #call regex library
import sys
import numpy as np
from scipy import constants as const
import subprocess
import readline, glob
import random
import math
import matplotlib.pyplot as plt

radius=20.0
peric=15.0

c=3.086e18
G=const.G
pc=const.parsec
Myr=const.mega*const.year
Msun=1.98855*10**30

nam_bin1_1=-99
nam_bin1_2=-99
nam_bin2_1=-99
nam_bin2_2=-99
m1=-99
m2=-99
m3=-99


######## tab completing input infos ########
def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]
readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

######## parabolic merging parameters calculator ########
def parabola(r, q):
    cos_theta=2*q/r -1
    sin_theta=(1-cos_theta**2)**0.5
    x=r*cos_theta
    y=r*sin_theta
    #v_x=0.25*sin_theta/np.sqrt(q)#ATTENZIONE!!!!!!! hai messo 1/4 di v parab!
    #v_y=-0.25*(cos_theta+1)/np.sqrt(q)#ATTENZIONE!!!!!!! hai messo 1/4 di v parab!
    v_x=sin_theta/np.sqrt(q)
    v_y=-(cos_theta+1)/np.sqrt(q)
    return x, y, v_x, v_y

#initial position shift calculator in order to have SCs on a filament, almost evenly spaced but with some degree of randomness (custom)
def filament(num, length):
    randomness_percent=1./10# 10% of randomness in position
    ic=np.ndarray(shape=(num, 3))
    side=length/float(np.sqrt(3))
    param=np.linspace(0.0, side, num=num)
    ic[:, 0] = param+np.random.random_sample((num,))*(side*randomness_percent)
    ic[:, 1] = param+np.random.random_sample((num,))*(side*randomness_percent)
    ic[:, 2] = param+np.random.random_sample((num,))*(side*randomness_percent)    
    return ic 

#initial position shift calculator in order to have SCs on a sphere on the x-y-z plane
def sphere(num, radius):
    ic=np.ndarray(shape=(num, 3))
    theta = (np.random.random_sample((num,))*2*math.pi)-math.pi
    cosphi = (np.random.random_sample((num,))*2)-1
    u = np.random.random_sample((num,))
    phi = np.arccos(cosphi)
    r = radius #*pow(u, float(1./3) ) # if one wants full sphere distrib
    ic[:, 0] = r * np.sin(phi) * np.cos( theta )
    ic[:, 1] = r * np.sin(phi) * np.sin( theta )
    ic[:, 2] = r * np.cos(phi)
    return ic 
    

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
                            
######## regular expression search ########
massa=re.compile('^  m  =\s+(\S+)')
mscal=re.compile('^  mass\_scale     =\s+(\S+)')
tscal=re.compile('^  time\_scale     =\s+(\S+)')
N_tot=re.compile('^  N =\s+(\S+)') #look at number of particles!
num_nodes=re.compile('-n\s(\d+)') #look at number of nodes(=centers of mass)#initial number of nodes
ind=re.compile('^  i =\s+(\S+)')
name=re.compile('^  name =\s+\((\d+),(\d+)\)')
name_sing=re.compile('^  name =\s+(\S+)')
name_merg=re.compile('^  name =\s+(\d+)\+(\d+)')
vel=re.compile('^  v  =\s+(\S+)\s(\S+)\s(\S+)')
rad=re.compile('^  r  =\s+(\S+)\s(\S+)\s(\S+)')
fine_head=re.compile('^\)Star')
com_t=re.compile('^  com_time =\s+(\S+)')
com_v=re.compile('^  com_vel =\s+(\S+)\s(\S+)\s(\S+)')
com_r=re.compile('^  com_pos =\s+(\S+)\s(\S+)\s(\S+)')
total_energy=re.compile('^  total_energy =\s+(\S+)')
begin=re.compile('^\(Particle')
log_fin=re.compile('\)Log')
log_in=re.compile('\(Log')
sys_time=re.compile('^  system_time  =\s+(\S+)')
time=re.compile('^  t  =\s+(\S+)')
dyn_in=re.compile('^  a  =\s+(\S+)\s(\S+)\s(\S+)')
dyn_fin=re.compile('\)Dynamics')
triple1=re.compile('^  name =\s+\((\d+),\((\d+),(\d+)\)\)')
triple2=re.compile('^  name =\s+\(\((\d+),(\d+)\),(\d+)\)')
triple3=re.compile('^  name =\s+\((\d+),(\d+),(\d+)\)')


######## I/O files request and opening ######## 
line=0
#num_sc_to_merge=int(raw_input('How many SCs to merge? '))
num_sc_to_merge=int(sys.argv[1])

f_in_name_=["" for a in range(num_sc_to_merge)]
f_in=[]
filend_=np.zeros(num_sc_to_merge)

#for numero in range(num_sc_to_merge):
    #f_in_name_[numero] = str(raw_input('Enter file name: '))
for numero in range(num_sc_to_merge):
    f_in_name_[numero] = str(sys.argv[2+numero])


#f_out_name=raw_input('Enter OUTPUT file:')
f_out_name=str(sys.argv[2+num_sc_to_merge])


for numero in range(num_sc_to_merge):
    if f_in_name_[numero].endswith(".gz"):
        f_in.append(gzip.open(f_in_name_[numero], "r"))
    else:
        f_in.append(open(f_in_name_[numero], "r"))
    filend_[numero]=file_len(f_in_name_[numero])

if f_out_name.endswith(".gz"):
    f_out=gzip.open(f_out_name, "w")
else:
    f_out=open(f_out_name, "w")


####### Find, store and reset some general system data: mscale, tscale, N, inital number of nodes, actual number of nodes... ####### 
m_=np.arange(num_sc_to_merge, dtype='float')
nodes_=np.arange(num_sc_to_merge)
nodes_actually_=np.arange(num_sc_to_merge)
msca_=np.arange(num_sc_to_merge, dtype='float')
tsca_=np.arange(num_sc_to_merge, dtype='float')
N_=np.arange(num_sc_to_merge)

for numero in range(num_sc_to_merge):
    line=0
    for s in f_in[numero]:
        line=line+1
        m_root=massa.search(s)
        msca=mscal.search(s)
        tsca=tscal.search(s)
        N=N_tot.search(s)
        node=num_nodes.search(s)
        if ((N != None) and ((line == 2)or(line == 3))): # total num of particles: line 2 or 3, because if the file is an t!=0 file, there will be the 'name=root' extra line
            N_[numero]=eval(N.group(1))

        if (node != None) and ((line == 5)or(line == 6)): # number of nodes: line 5 or 6 as above
            nodes_[numero]=eval(node.group(1))
            
        if (m_root != None): 
            m_[numero]=eval(m_root.group(1))
            
        if (msca != None):
            msca_[numero]=eval(msca.group(1))
            
        if(tsca != None):
            tsca_[numero]=eval(tsca.group(1))
            break

#seek(0,0) comes back to the beginning of clusters
    f_in[numero].seek(0,0)
    
#search for the actual number of nodes = N (number of particles) - ((numbr of lines N!= 1)-1)
#e.g. for evey binary system (N = 2) it will subtract 1 to the number of particles to have the actual number of nodes
for numero in range(num_sc_to_merge):
    line=0
    count_bin_tr=0
    for s in f_in[numero]:
        N=N_tot.search(s)
        line=line+1
        if ((N != None) and ((line != 2)and(line != 3))):
            if eval(N.group(1)) != 1:
                count_bin_tr+=eval(N.group(1))-1
    nodes_actually_[numero] = N_[numero] - count_bin_tr

#seek(0,0) comes back to the beginning of clusters
    f_in[numero].seek(0,0)
#set new mscale (it is called new_mtot, but it is actually new_mscale), tscale, nodes etc values
print 'initial number of particles: ', N_
print 'initial number of nodes: ', nodes_
print 'actual number of nodes: ', nodes_actually_
#new_mtot = 1.0/sum(1./msca_)#old definition of new_mtot, it works if m_root is ~1, i.e. if you are merging std input files
new_mtot = 1.0/sum(m_/msca_)#new definition of new_mtot, more general, it works even if you are merging std input files and evolved snapshots
new_N = sum(N_)
new_nodes = sum(nodes_)
new_nodes_actually =sum(nodes_actually_)
new_tsca=tsca_[0]*np.sqrt(msca_[0]/new_mtot)

#power=len(str(new_nodes)) # needed for binary names # this line is correct if you are merging two std input files, otherwise it is better to consider the actual number of nodes, which can in principle be different from the initial one in an file already evolved by kira
power=len(str(new_nodes_actually)) # needed for binary names

########    get shift in position and velocity to merge ########
pos_shift=np.ndarray(shape=(num_sc_to_merge,3))
v_rel=np.ndarray(shape=(num_sc_to_merge,3))

v_rel[:][:]=0

#configuration = str(raw_input('Select initial spatial configuration: filament or sphere?'))
configuration = str(sys.argv[2+num_sc_to_merge+1])
#extent = int(raw_input('Select configuration extent (in pc):'))
extent = int(sys.argv[2+num_sc_to_merge+2])
if configuration == 'filament':
    pos_shift=filament(num_sc_to_merge, extent)
if configuration == 'sphere':
    pos_shift=sphere(num_sc_to_merge, extent)
if (configuration != 'sphere') and (configuration != 'filament'):
    print "Configuration not available!"


#################################################################
#                         FILE 1                                #
#################################################################


#put mscale to -1 (to avoid problem with   m  =  1.0000 after system time)
mscale=-1.0
line=0
log_switch_on=0
log_switch_off=0
log_trig=0
dyn_switch=0
root_info=0
#modify and rewrite file1
for s in f_in[0]:
    
    line=line+1
    r=rad.search(s)
    cm_r=com_r.search(s)
    cm_v=com_v.search(s)
    cm_t=com_t.search(s)
    msca=mscal.search(s)
    masse=massa.search(s)
    tsca=tscal.search(s)
    nam=name.search(s)    
    N=N_tot.search(s)
    i=ind.search(s)
    v=vel.search(s)
    beg=begin.search(s)
    f_log=log_fin.search(s)
    i_log=log_in.search(s)
    t_e=total_energy.search(s)
    s_t=sys_time.search(s)
    nam_merg=name_merg.search(s)
    t=time.search(s)
    i_dyn=dyn_in.search(s)
    f_dyn=dyn_fin.search(s)
    nam_sing=name_sing.search(s)
    tr1=triple1.search(s)
    tr2=triple2.search(s)
    tr3=triple3.search(s)

#eliminates initial info on the root particle: anyway kira calculates them correctly at the first run
    if mscale <0.0 and (masse !=None or v!=None or r!=None or cm_r!=None or cm_t!=None or cm_v!=None):
        root_info=1
    else: root_info=0

#eliminate some lines in Dynamics, eventually produced by kira 
    if (i_dyn != None):
        dyn_switch=1

    if (f_dyn != None):
        dyn_switch=0

#set system time to 0
    if (s_t != None):
        if (eval(s_t.group(1)) != 0):
            f_out.write('  system_time  =  0'+'\n')
            
#add 'i=0' just before the total number of particles to make kira compile it
    if ((beg != None) and (line == 1)):
        f_out.write(s)
        f_out.write('  i = 0'+'\n')
        
#modify only the first occurance of N (line 2 or 3: see above)  
    if ((N != None) and (line == 2 or line ==3)): 
        f_out.write('  N = '+str(new_N)+'\n')
        
#delete the content in first (Log-Log) using some switches 
    if log_trig==1:
        log_switch_on=1
    if i_log != None:
        log_trig=1
    if f_log!=None:
        log_switch_off=1

#read old masscale and set new one (same for timescale)
    if (msca != None):
        mscale=eval(msca.group(1))
        f_out.write('  mass_scale     =  '+str(new_mtot)+'\n')
    if (tsca != None):
        f_out.write('  time_scale     =  '+str(new_tsca)+'\n')

#update 'name' and index in case of binary
    if (nam != None):
        nam_bin1_1=eval(nam.group(1))
        nam_bin1_2=eval(nam.group(2))
        if nam_bin1_1 <= nodes_[0]:
            indice_bin1=nam_bin1_1
        else:
            indice_bin1=10**power+(nam_bin1_1 % 10**(len(str(nodes_[0]))))
        if nam_bin1_2 <= nodes_[0]:
            indice_bin2=nam_bin1_2
        else:
            indice_bin2=10**power+(nam_bin1_2 % 10**(len(str(nodes_[0]))))
        f_out.write('  name = ('+str(indice_bin1)+','+str(indice_bin2)+')'+'\n')
        
            
#update 'name' for merged stars, only in the case they come from binary system (else just copy old names)
    if (nam_merg !=None):
        part1 = eval(nam_merg.group(1))
        part2 = eval(nam_merg.group(2))
        if part1 <= nodes_[0]:
                progen_1=part1
        else:
                progen_1= 10**power+(part1 % 10**(len(str(nodes_[0]))))
        if part2 <= nodes_[0]:
                progen_2=part2
        else:
                progen_2= 10**power+(part2 % 10**(len(str(nodes_[0]))))

        f_out.write('  name = '+str(progen_1)+'+'+str(progen_2)+'\n')


#update triples (assumed the possible configurations are the following: ((bin1,bin2),star), (star,(bin1,bin2)), ((bin1,star),bin2), (bin1,(star,bin2)), (star1,star2,star3) or obvious permutations) and anyway a warning message is given

    if (tr1 != None) or (tr2 != None) or (tr3 != None) :
        print "Warning! Triple system found! Check out!"
        m1=eval(tr1.group(1))
        m2=eval(tr2.group(2))
        m3=eval(tr3.group(3))
        if m1 <= nodes_[0]:
            new_m1=m1
        else:
            new_m1=10**power+(m1 % 10**(len(str(nodes_[0]))))
        if m2 <= nodes_[0]:
            new_m2=m2
        else:
            new_m2=10**power+(m2 % 10**(len(str(nodes_[0]))))
        if m3 <= nodes_[0]:
            new_m3=m3
        else:
            new_m3=10**power+(m3 % 10**(len(str(nodes_[0]))))

        if tr1 != None:
            f_out.write('  name = ('+str(new_m1)+',('+str(new_m2)+','+str(new_m3)+'))'+'\n')
        if tr2 != None:
            f_out.write('  name = (('+str(new_m1)+','+str(new_m2)+'),'+str(new_m3)+')'+'\n')
        if tr3 != None:
            f_out.write('  name = ('+str(new_m1)+','+str(new_m2)+','+str(new_m3)+')'+'\n')
                
#update index 'i' in case of binary/triple companion, else just copy the old one  
    if (i != None):
        index=eval(i.group(1))
        if (index <= nodes_[0]): #index i depends first on the numbr of nodes
            indice=index
        else:
            indice = 10**power+(index % 10**(len(str(nodes_[0]))))
        f_out.write('  i = ' +str(indice)+'\n')
                
#update single stars masses according to new_mtot
    if((masse != None) and (mscale>0.0)):
        m=(eval(masse.group(1))/mscale)*new_mtot
        f_out.write('  m  =  '+str(m)+'\n')
        
#update single stars velcoities according to new_tsca
    if ((v != None) and (mscale>0.0)): #we are ignoring the first v, referring to he system cm velocity 
        v_x=(eval(v.group(1))*tsca_[0])/new_tsca
        v_y=(eval(v.group(2))*tsca_[0])/new_tsca
        v_z=(eval(v.group(3))*tsca_[0])/new_tsca
        f_out.write('  v  =  '+str(v_x)+' '+str(v_y)+' '+str(v_z)+'\n')
        
#check the values of Log switches to cut off the content 
    log_status=log_switch_on-log_switch_off
    if log_status == 0 and log_switch_off==1:
        log_switch_on=0
        log_switch_off=0
        log_trig=0
            
#writes the lines that must remain UNCHANGED in the new file
    if(((masse == None) and (root_info == 0) and (msca == None) and (t == None) and (tsca == None) and (tr1 == None) and (tr2 == None) and (tr3 == None) and (dyn_switch == 0) and (s_t == None) and (i == None) and (N == None) and (nam == None) and (nam_merg == None) and (v == None) and (beg == None) and (t_e == None) and (nam_sing == None) and (log_status == 0)) or (t_e ==None and log_status == 0 and line!=2 and line!=3 and line!=1 and nam == None and (tr1 == None) and (tr2 == None) and (tr3 == None) and s_t == None and dyn_switch == 0 and root_info == 0 and t == None and (mscale<0.0 or N!= None or beg!=None))):
        if(line < filend_[0]):
            f_out.write(s) #copies line to output file 


#################################################################
#                      FILE 2, 3, 4, 5 ...                      #
#################################################################

for numero in range(1,num_sc_to_merge):
    line=0
    dyn_switch=0
    mscale=-1.0
    indice_bin1=-99
    indice_bin2=-99
    indice=-999
    new_m1=-99
    new_m2=-99
    new_m3=-99
    
#this delete file 2's header 
    for s in f_in[numero]:
        line=line+1
        fin=fine_head.search(s)
        msca=mscal.search(s)
        if(msca != None):
            mscale=eval(msca.group(1))
        if(fin!=None):
            break
        
#now attach file 2, 3, 4... (some stuff identical as before)
    for s in f_in[numero]:
        line=line+1 
        masse=massa.search(s)
        i=ind.search(s)
        nam=name.search(s)
        nam_merg=name_merg.search(s)
        nam_sing=name_sing.search(s)
        v=vel.search(s)
        r=rad.search(s)
        fin=fine_head.search(s)
        s_t=sys_time.search(s)
        t=time.search(s)
        i_dyn=dyn_in.search(s)
        f_dyn=dyn_fin.search(s)
        tr1=triple1.search(s)
        tr2=triple2.search(s)
        tr3=triple3.search(s)

#eliminate some lines in Dynamics, eventually produced by kira 
        if (i_dyn != None):
            dyn_switch=1

        if (f_dyn != None):
            dyn_switch=0
#set system time to 0
        if (s_t != None):
            if (eval(s_t.group(1)) != 0):
                f_out.write('  system_time  =  0'+'\n')

#update name in the binaries
        if (nam != None):
            nam_bin2_1=eval(nam.group(1))
            nam_bin2_2=eval(nam.group(2))
            if nam_bin2_1 <= nodes_[numero]:
                indice_bin1=nam_bin2_1+sum(nodes_[:numero])
            else:
                indice_bin1=10**power+(nam_bin2_1 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
            if nam_bin2_2 <= nodes_[numero]:
                indice_bin2=nam_bin1_2+sum(nodes_[:numero])
            else:
                indice_bin2=10**power+(nam_bin2_2 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
            f_out.write('  name = ('+str(indice_bin1)+','+str(indice_bin2)+')'+'\n')
 
#update 'name' for merged stars, only in the case the come from binary system (else just copy old names)

        if (nam_merg !=None):
            part1 = eval(nam_merg.group(1))
            part2 = eval(nam_merg.group(2))
            if part1 <= nodes_[numero]:
                progen_1=part1+sum(nodes_[:numero])
            else:
                progen_1= 10**power+(part1 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
            if part2 <= nodes_[numero]:
                progen_2=part2+sum(nodes_[:numero])
            else:
                progen_2= 10**power+(part2 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])

            f_out.write('  name = '+str(progen_1)+'+'+str(progen_2)+'\n')

#update triples (assumed the possible configurations are the following: ((bin1,bin2),star), (star,(bin1,bin2)), ((bin1,star),bin2), (bin1,(star,bin2)), (star1,star2,star3) or obvious permutations) and anyway a warning message is given

            if (tr1 != None) or (tr2 != None) or (tr3 != None) :
                print "Warning! Triple system found! Check out!"
                m1=eval(tr1.group(1))
                m2=eval(tr2.group(2))
                m3=eval(tr3.group(3))
                if m1 <= nodes_[numero]:
                    new_m1=m1+sum(nodes_[:numero])
                else:
                    new_m1=10**power+(m1 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
                if m2 <= nodes_[numero]:
                    new_m2=m2+sum(nodes_[:numero])
                else:
                    new_m2=10**power+(m2 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
                if m3 <= nodes_[numero]:
                    new_m3=m3+sum(nodes_[:numero])
                else:
                    new_m3=10**power+(m3 % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])

                if tr1 != None:
                    f_out.write('  name = ('+str(new_m1)+',('+str(new_m2)+','+str(new_m3)+'))'+'\n')
                if tr2 != None:
                    f_out.write('  name = (('+str(new_m1)+','+str(new_m2)+'),'+str(new_m3)+')'+'\n')
                if tr3 != None:
                    f_out.write('  name = ('+str(new_m1)+','+str(new_m2)+','+str(new_m3)+')'+'\n')           

#update single particles index i, checking if it is a stand alone particle or binary partner

        if (i != None):
            index=eval(i.group(1))
            if (index <= nodes_[numero]): #index i depends first on the numbr of nodes
                    indice=index+sum(nodes_[:numero])
            else:
                indice = 10**power+(index % 10**(len(str(nodes_[numero]))))+sum(nodes_[:numero])
            f_out.write('  i = ' +str(indice)+'\n')
            
            
#update single masses
        if((masse != None) and (mscale>0.0)):
            m=(eval(masse.group(1))/mscale)*new_mtot
            f_out.write('  m  =  '+str(m)+'\n')

#convert velocity to new units(timescale) AND ONLY for NOT binary members add a shift in velocity, to have parabolic orbit (NB: it skips the velocity of root particle - very first occurence of v-)
        if ((v != None) and (mscale>0.0)):
            v_x=(eval(v.group(1))*tsca_[numero])/new_tsca
            v_y=(eval(v.group(2))*tsca_[numero])/new_tsca
            v_z=(eval(v.group(3))*tsca_[numero])/new_tsca
            if (indice != indice_bin1 and indice != indice_bin2 and indice != new_m1 and indice != new_m2 and indice != new_m3):#this is the line that check we are not in binary
                v_x_s=v_x+v_rel[numero][0]
                v_y_s=v_y+v_rel[numero][1]
                v_z_s=v_z+v_rel[numero][2]
                f_out.write('  v  =  '+str(v_x_s)+' '+str(v_y_s)+' '+str(v_z_s)+'\n')
            else:
                f_out.write('  v  =  '+str(v_x)+' '+str(v_y)+' '+str(v_z)+'\n')
      
#gives a shift in position (ONLY for NOT binary members-meaning it doesn't change internal coordinates, just CM of the binary-)
        if (indice != indice_bin1 and indice != indice_bin2 and indice != new_m1 and indice != new_m2 and indice != new_m3):
            if ((r!=None) and (mscale>0.0)):
                r_x=eval(r.group(1))+pos_shift[numero][0]
                r_y=eval(r.group(2))+pos_shift[numero][1]
                r_z=eval(r.group(3))+pos_shift[numero][2]
                f_out.write('  r  =  '+str(r_x)+' '+str(r_y)+' '+str(r_z)+'\n')
        else:
            if ((r!=None) and (mscale>0.0)):
                f_out.write(s)


        
#copies all the rest unchanged
        if(((masse == None) and (msca == None) and (dyn_switch == 0) and (tr1 == None) and (tr2 == None) and (tr3 == None) and (tsca == None) and (t == None) and (s_t == None) and (i == None) and (nam_merg == None) and (nam_sing == None)  and (nam == None) and (v == None) and (r == None)) or (nam == None and mscale < 0.0)):
            if line < filend_[numero]:
                f_out.write(s)
            if (numero == (num_sc_to_merge-1))*(line == filend_[numero]):
                f_out.write(s)
        if (fin!=None):
            continue
    
#close files
for numero in range(0,num_sc_to_merge):
    f_out.close()
    f_in[numero].close()
    


  
 
