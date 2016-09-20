"""
Library of function declearations used in plotting verificaiton results
"""
import scipy as sp
import matplotlib
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import cm
from scipy import linalg
from shutil import copyfile
import os 
import sys

class ThermoData(object):
    Temperature=False;
    Density=False;
    Enthalpy=False;
    StaticEnergy=False;
    Entropy=False;
    Pressure=False;
    SoundSpeed2=False;
    dPdrho_e=False;
    dPde_rho=False;
    dTdrho_e=False;
    dTde_rho=False;
    Cp=False;
    Mu=False;
    Kt=False;
    
    def load_data(self, filename):
        self.data = sp.genfromtxt(filename, delimiter=', ')
        self.Temperature  = self.data[:,0];
        self.Density      = self.data[:,1];
        self.Enthalpy     = self.data[:,2];
        self.StaticEnergy = self.data[:,3];
        self.Entropy      = self.data[:,4];
        self.Pressure     = self.data[:,5];
        self.SoundSpeed2  = self.data[:,6];
        self.dPdrho_e     = self.data[:,7];
        self.dPde_rho     = self.data[:,8];
        self.dTdrho_e     = self.data[:,9];
        self.dTde_rho     = self.data[:,10];
        self.Cp           = self.data[:,11];
        self.Mu           = self.data[:,12];
        #self.dmudrho_T    = self.data[:,13];
        #self.dmudT_rho    = self.data[:,14];
        self.Kt           = self.data[:,13];
        #self.dktdrho_T    = self.data[:,16];
        #self.dktdT_rho    = self.data[:,17];
        return 
        
    def load_FP_output(self, filename):
        self.data = sp.genfromtxt(filename, skip_header=2)
        self.Density      = self.data[:,0];
        self.Pressure     = self.data[:,1];
        self.SoundSpeed2  = self.data[:,2]**2;
        self.Cp           = self.data[:,3];
        self.Entropy      = self.data[:,4];
        self.Mu           = self.data[:,5];
        self.Kt           = self.data[:,6];
        self.dPdrho_e     = self.data[:,7];
        self.dPde_rho     = self.data[:,8];
        self.dTdrho_e     = self.data[:,9];
        self.dTde_rho     = self.data[:,10];
        self.Temperature  = self.data[:,11];
        self.StaticEnergy = self.data[:,12];
        self.Enthalpy     = self.data[:,13];
        return;
    def load_rand_output(self, filename):
        self.data = sp.genfromtxt(filename, skip_header=1)
        self.Density      = self.data[:,0];
        self.Pressure     = self.data[:,1];
        self.SoundSpeed2  = self.data[:,2]**2;
        self.Cp           = self.data[:,3];
        self.Entropy      = self.data[:,4];
        self.Mu           = self.data[:,5];
        self.Kt           = self.data[:,6];
        self.dPdrho_e     = self.data[:,7];
        self.dPde_rho     = self.data[:,8];
        self.dTdrho_e     = self.data[:,9];
        self.dTde_rho     = self.data[:,10];
        self.Temperature  = self.data[:,11];
        self.StaticEnergy = self.data[:,12];
        self.Enthalpy     = self.data[:,13];
        
    def plot_2D(self, axis_x, axis_y):
        plt.plot(getattr(self, axis_x), getattr(self,axis_y),'ro',alpha=0.1);
        return
        

class RandomSamples(ThermoData):
    
    def __init__(self,filename):
        print 'Loading random verification data'
        self.load_rand_output(filename);
        print 'DONE Loading random verification data'
        
        #Prepare the input files
        print 'Preparing input files...'
        sp.savetxt('rhoe_in.dat',sp.column_stack((self.Density,self.StaticEnergy)), delimiter='\t')
        sp.savetxt('PT_in.dat'  ,sp.column_stack((self.Pressure,self.Temperature)), delimiter='\t')
        sp.savetxt("Prho_in.dat",sp.column_stack((self.Pressure,self.Density)), delimiter='\t')
        sp.savetxt("rhoT_in.dat",sp.column_stack((self.Density,self.Temperature)), delimiter='\t')
        sp.savetxt("Ps_in.dat"  ,sp.column_stack((self.Pressure,self.Entropy)), delimiter='\t')
        sp.savetxt("hs_in.dat"  ,sp.column_stack((self.Enthalpy,self.Entropy)), delimiter='\t')        
        print 'DONE Preparing input files...'
        return


class RestartSamples(object):
    Density=None;
    Pressure=None;
    Temperature=None;
    restarts = None;
    
    def __init__(self,filename, dim='2D'):
        restarts = sp.genfromtxt(filename, skip_header=2)
        if dim == '2D':
            self.Density      = restarts[:,3];
            self.Pressure     = restarts[:,9];
            self.Temperature  = restarts[:,10];        
        elif dim == '3D':
            self.Density      = restarts[:,4];
            self.Pressure     = restarts[:,9];
            self.Temperature  = restarts[:,10];        
        self.restarts=restarts;
            
        print 'Restart points: ', len(self.Density)
        return;
        
    def plot_restarts(self, scatter_x ='Density', scatter_y='Pressure'):        
        plt.plot(getattr(self,scatter_x), getattr(self,scatter_y),'.',c='#B2FF66',label='PR_gas Restart');
        return;


        
class SciPy_InterpolatedData(ThermoData):
    median_ERR = 1;
    

    def __init__(self, which_case, LUT, RandomSamples, interp_type):
        print 'SciPy Interpolating ', which_case
        
        select = {\
        "rhoe":('Density','StaticEnergy'),\
        "PT":('Pressure','Temperature'),\
        "Prho":('Pressure','Density'),\
        "rhoT":('Density','Temperature'),\
        "Ps":('Pressure','Entropy'),\
        "hs":('Enthalpy','Entropy')\
        }
        
        thermo1, thermo2, = select[which_case]
        x =getattr(LUT,thermo1)
        y =getattr(LUT,thermo2)
        samples_x = getattr(RandomSamples,thermo1)
        samples_y = getattr(RandomSamples,thermo2)
        setattr(self,thermo1, samples_x)
        setattr(self,thermo2, samples_y)
        
        #variables = sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
        #'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho',\
        #'dTdrho_e','dTde_rho','Cp'])#,'Mu','Kt']);
        #variables = sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
        #'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho','Cp'])#,'Mu','Kt']);
        
        for var in variables[sp.where((variables!=thermo1) * (variables!=thermo2))]:
            z = getattr(LUT,var)            
            interp_func = sp.interpolate.griddata((x,y),z,sp.column_stack((samples_x,samples_y)),\
            method=interp_type) 
            nan_index = sp.where(sp.isnan(interp_func))
            interp_func[nan_index]= sp.interpolate.griddata((x,y),z,\
            sp.column_stack((samples_x[nan_index],samples_y[nan_index])),\
            method='nearest') 
            setattr(self,var,interp_func)
            
        return  
          
        
class SU2_InterpolatedData(ThermoData):
    name='';
    data=False;
    median_ERR = 1;
    
    def __init__(self, filename):
        self.load_data(filename);

        
class PRGrid(ThermoData):
    
    meshfile=False;
    P_dim=0;
    D_dim=0;
        
    def __init__(self, filename):
        self.load_data(filename);
        self.P_dim = sp.where((self.Density-self.Density[0])!=0)[0][0]
        self.D_dim = len(self.Density)/self.P_dim
        print 'P-RHO dimensions: %i by %i'%(self.P_dim, self.D_dim);

    def plot_mesh(self,thermo_x,thermo_y):
        x = getattr(self,thermo_x)
        y = getattr(self,thermo_y)
        
        for i in range(self.P_dim):
            plt.plot(x[i:][::self.P_dim], y[i:][::self.P_dim], 'c-', alpha=0.3);
        for i in range(self.D_dim):
            plt.plot(x[i*self.P_dim:(i+1)*self.P_dim], y[i*self.P_dim:(i+1)*self.P_dim], 'c-', alpha=0.3);

class PRGridSkewed(ThermoData):
    
    meshfile=False;
    P_dim=0;
    D_dim=0;
        
    def __init__(self, filename):
        self.load_FP_output(filename)
        self.P_dim = sp.where((self.Density-self.Density[0])!=0)[0][0]
        self.D_dim = len(self.Density)/self.P_dim
        print 'P-RHO dimensions: %i by %i'%(self.P_dim, self.D_dim);

    def plot_mesh(self,thermo_x,thermo_y):
        x = getattr(self,thermo_x)
        y = getattr(self,thermo_y)
        
        for i in range(self.P_dim):
            plt.plot(x[i:][::self.P_dim], y[i:][::self.P_dim], 'k-', alpha=0.5);
        for i in range(self.D_dim):
            plt.plot(x[i*self.P_dim:(i+1)*self.P_dim], y[i*self.P_dim:(i+1)*self.P_dim], 'k-', alpha=0.5);    
        plt.grid(which='both')
           
    

class RefinementLevel(object):
   cases=False;
   filename=False;
   #variables=sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
    #    'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho',\
     #   'dTdrho_e','dTde_rho','Cp','Mu','Kt']);
   #variables=sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
    #   'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho','Cp','Mu']);
   variables=sp.array(['Temperature','Density','Enthalpy','StaticEnergy',\
       'Entropy','Pressure','SoundSpeed2','dPdrho_e','dPde_rho','Cp']);
       
   select = {\
        "rhoe":('Density','StaticEnergy'),\
        "PT":('Pressure','Temperature'),\
        "Prho":('Pressure','Density'),\
        "rhoT":('Density','Temperature'),\
        "Ps":('Pressure','Entropy'),\
        "hs":('Enthalpy','Entropy')\
        }
   LUT={};
   SU2={};
   time=None;
        
   RandomSamples=False;
   SciPy={}; 
            
   def __init__(self, filename, cases):        
       self.cases =cases;        
       self.filename=filename;
       self.time = {"rhoe":1000,"PT":1000,"Prho":1000,\
        "rhoT":1000,"Ps":1000,"hs":1000};
   def load_random_samples(self,RandomSamples):
       self.RandomSamples = RandomSamples
   
   def load_mesh(self):
       self.LUT = PRGrid('mesh.dat')
       
   def load_results_SU2(self):
       self.SU2={i: SU2_InterpolatedData(str(i+'_out.dat')) for i in self.cases};            
            
   def load_results_SciPy(self,interp_kind='linear'):
       self.SciPy={i:SciPy_InterpolatedData(i,self.LUT, self.RandomSamples, \
       interp_kind ) for i in self.cases};
     
   def load_time_perf(self,filename):
       with open(filename, 'r') as f:
           lines = f.readlines()
       for j in self.time.keys():
           for i in range(len(lines)):
               if (lines[i].find(j) != -1):
                   print lines[i]
                   self.time[j] = float(lines[i].strip().split()[-1])
               
      
   def get_REL_ERR_SU2(self,which_case):
        i=0;
        thermo1 = self.select[which_case][0]
        thermo2 = self.select[which_case][1]
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SU2[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR)/i
        setattr(self.SU2[which_case],"median_ERR",sp.median(self.REL_ERR));
        return
    
   def get_REL_ERR_SciPy(self,which_case):
        i=0;
        thermo1 = self.select[which_case][0]
        thermo2 = self.select[which_case][1]
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SciPy[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR)/i
        setattr(self.SciPy[which_case],"median_ERR",sp.median(self.REL_ERR));
        return
       
       
   def plot_REL_ERR_SU2(self,which_case):
        i=0;
        thermo1 = self.select[which_case][0]
        thermo2 = self.select[which_case][1]
        self.get_REL_ERR_SU2(which_case)
        
        print 'Median error SU2', sp.median(self.REL_ERR)
        print 'Mean error SU2', sp.mean(self.REL_ERR)
        print 'Max error SU2', max(self.REL_ERR)
        print 'Min error SU2', min(self.REL_ERR)
        x = getattr(self.SU2[which_case],thermo1)
        y = getattr(self.SU2[which_case],thermo2)
        #trusted_values = sp.where(self.REL_ERR>0<0.9*max(self.REL_ERR))
        #self.REL_ERR = self.REL_ERR[trusted_values]
        #x = x[trusted_values]
        #y = y[trusted_values]
        scat=plt.scatter(x,y,c=self.REL_ERR)                
        plt.grid(which='both')
        scat.set_array(self.REL_ERR)        
        plt.colorbar(scat)
        plt.xlim((min(x)*0.95,max(x)*1.05));
        plt.ylim((min(y)*0.95,max(y)*1.05));
        print 'x argmax %i , x_val: %f ' %(sp.argmax(self.REL_ERR),x[sp.argmax(self.REL_ERR)])
        print 'y argmax %i , y_val: %f ' %(sp.argmax(self.REL_ERR),y[sp.argmax(self.REL_ERR)])
        return;
       
   
   def plot_hist_compare(self,which_case):
        plt.ylabel('Percentage of points')
        plt.xlabel('Percentage RMS relative error')
        
        def yto_percent(y, x):
            s = str(sp.around((y/(len(self.REL_ERR)*1.0)*100),2))
            if matplotlib.rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'     

        def xto_percent(y, x):
            s = str(y*100)
            if matplotlib.rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                    return s + '%' 
        
        thermo1, thermo2, = self.select[which_case]
        #Plot the SU2 error
        i=0;
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SU2[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR)/i
        plt.hist(self.REL_ERR, bins=25, color='k', alpha=0.3, label='SU2')
        print 'Error max SU2', max(self.REL_ERR)
        setattr(self.SU2[which_case],"median_ERR",sp.median(self.REL_ERR));
        
        #Plot the SciPy error
        i =0;
        self.REL_ERR = 0;
        for v in self.variables[sp.where\
        ((self.variables!=thermo1) * (self.variables!=thermo2))]:
            i=i+1;
            self.REL_ERR = self.REL_ERR + \
            ((getattr(self.SciPy[which_case],v)-getattr(self.RandomSamples,v))/\
            (getattr(self.RandomSamples,v)))**2;
        self.REL_ERR = sp.sqrt(self.REL_ERR)/i
        
        plt.hist(self.REL_ERR, bins=25, color='c', alpha=0.5, label='SciPy')
        print 'Error max SciPy', max(self.REL_ERR)
        setattr(self.SciPy[which_case],"median_ERR",sp.median(self.REL_ERR));

        
        formatter_y = FuncFormatter(yto_percent)
        formatter_x = FuncFormatter(xto_percent)
        plt.gca().yaxis.set_major_formatter(formatter_y)
        plt.gca().xaxis.set_major_formatter(formatter_x)
        plt.grid(which='both')
        plt.legend()

       
        return       
        
def plot_time_perf(RefinementLevels, FP_time, points_count):
    for i in RefinementLevels[0].time.keys():
            x =[];
            y =[];
            print "Analyzing computation time of: ", i ;
            for r in RefinementLevels:                
                if r!= RefinementLevels[0]:
                    x.append(r.LUT.D_dim*r.LUT.P_dim)
                    y.append(r.time[i])
            x = sp.array(x)
            y = sp.array(y)            
            y = y[sp.argsort(x)]
            x = x[sp.argsort(x)]
            print y            
                                    
            plt.loglog(x,y/points_count, label='%s'%i, basex=10, basey=10) 
    plt.loglog(x,FP_time*sp.ones_like(x)/points_count,label='FluidProp',basex=10, basey=10)
    plt.grid(which='both')
    plt.xlabel('LuT Grid Nodes (N)')
    plt.ylabel('Time per point [sample points: %i]'%points_count)
    return
    
def plot_relative_perf(RefinementLevels, FP_time, points_count):
    for i in RefinementLevels[0].time.keys():
            x =[];
            y =[];
            print "Analyzing relative computation time of: ", i ;
            for r in RefinementLevels:  
                if r!= RefinementLevels[0]:
                    x.append(r.LUT.D_dim*r.LUT.P_dim)
                    y.append(r.time[i])
            print y
            x = sp.array(x)
            y = sp.array(y)            
            y = y[sp.argsort(x)]
            x = x[sp.argsort(x)]
                                    
            plt.semilogx(x,(FP_time/y), label='%s'%i)#, basex=10, basey=10) \
            #           subsy=sp.linspace(10**(-5), 10**(-2),20),\
                       #subsx=sp.linspace(10**(2), 10**(5),50))
    #plt.loglog(x,FP_time*sp.ones_like(x)/points_count,label='FluidProp',basex=10, basey=10)
    plt.grid(which='both')
    plt.xlabel('LuT Grid Nodes (N)')
    plt.ylabel('Speed-up factor wrt. FluidProp')
    return
    
        
def plot_median_errors(RefinementLevels):
        for i in RefinementLevels[0].cases:
            x =[];
            y =[];
            print "Analyzing median error on: ", i ;
            for r in RefinementLevels:                
                x.append(r.LUT.D_dim*r.LUT.P_dim)
                r.get_REL_ERR_SU2(i)
                y.append(r.SU2[i].median_ERR*100)
            
            x = sp.array(x)
            y = sp.array(y)            
            y = y[sp.argsort(x)]
            x = x[sp.argsort(x)]
                                    
            LHM = sp.ones((len(x),2))
            RHS = sp.ones((len(x),1))            
            LHM[:,1] = sp.log10(x)
            RHS[:,0] = sp.log10(y)

            sols = sp.linalg.lstsq(LHM,RHS)
            b = -sols[0][1]
            plt.loglog(x,y, label='%s, %s'%(i,r'$O(\frac{1}{N})^{%s}$'%str(sp.around(b,2))), basex=10, basey=10, \
                       subsy=sp.linspace(10**(-5), 10**(-2),20),\
                       subsx=sp.linspace(10**(2), 10**(5),50))
            
            #for r in RefinementLevels:                
               # x.append(r.LUT.D_dim*r.LUT.P_dim)
              #  r.get_REL_ERR_SciPy(i)
             #   y.append(r.SciPy[i].median_ERR*100)
            #plt.plot(x,y, label='SciPy: %s'%i)
        plt.grid(which='both')
        plt.xlabel('Grid Nodes (N)')
        plt.ylabel('Median relative error [%]')
        return;
       
              
if __name__ == '__main__':
    PR_files = os.listdir('TableLibrary')
    levels = len(PR_files);
    cases = ["rhoe","PT", "Prho","rhoT","Ps","hs"]
    Refinement_Levels = [RefinementLevel(PR_files[i], cases) for i in range(levels)]
    Random_Samples = RandomSamples('random.dat')
    for level in Refinement_Levels:
        output_files = ["rhoe_out.dat","PT_out.dat", "Prho_out.dat","rhoT_out.dat","Ps_out.dat","hs_out.dat"]
        for f in output_files:
            try:
                 os.remove(f)
            except OSError:
                 pass
             
        #Use the thermotable related to the desired refinement level
        copyfile('TableLibrary/'+level.filename, 'CO2.rgp')  
        print 'Mesh file:',  level.filename;
        #Recompile if necessary             
        os.system('make Debug/makefile');
        print 'Running...'
        #Run and save the log
        os.system('Debug/TableReader>log');
        print 'DONE'
        level.load_mesh();
        level.load_random_samples(Random_Samples)
        level.load_results_SU2();
        level.load_results_SciPy('linear');
        #level.LUT.plot_mesh('Density','StaticEnergy')
        plt.figure(figsize=(16,6))
        level.plot_hist_compare('rhoe')
        #level.plot_hist_compare('hs')
        #level.plot_hist_compare('Prho')


