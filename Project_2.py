import numpy as np
import matplotlib.pyplot as plt
import copy


#%% 1)


## We create the check function, which verifies that the random particle 
#  generated does not overlap with the rest of the particles



def check(x, y, particulas, diametro=1):   
    cond = all(
        [(np.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2) >= sigma) for p in particulas]
    )
    return cond

np.random.seed()

N = 1000            ## Number of particles
sigma = 1           ##  Diameter of the particles
L = np.sqrt(np.pi*N/(4*0.05))       ## Box' side length

## Particulas = [[x1, y1], [x2, y2], etc...]
particulas_i = []                 ## Array of particles
marker_size = np.pi * sigma**2  ## To plot the particles with a diameter equal to sigma


## We create a while loop that creates particles randomly distributed within the box
#  and we include the particles who do not overlap to our system and discard
#  those that overlap.

while len(particulas_i) < 1000:

    x = np.random.rand() * (L-sigma) - (L-sigma)/2
    y = np.random.rand() * (L-sigma) - (L-sigma)/2

    if check(x, y, particulas_i, diametro=sigma):
        particulas_i.append([x, y])

## We plot the particles to check that they, indeed, do not overlap

particulas_i = np.array(particulas_i)
plt.plot(particulas_i[:, 0], particulas_i[:, 1], "or", markersize=marker_size)
plt.show()





## We need to move all the particles, each one in a different random direction
#  we 

s = []

delta1 = 0.001            ## The amplitude of the random motion

particulas1 = copy.deepcopy(particulas_i)  ## We save the initial positions

dt = 1               ## Time jump, the unit of time is one Monte Carlo step

T = 100*dt          ## Total number of Monte Carlo steps


t = np.linspace(0, T, int(T/dt))        ## Our time

MSD1 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    
      while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*N)           ## Index to choose a random particle
    
        trial = copy.deepcopy(particulas1)
    
        x = particulas1[d][0] + (np.random.rand()-0.5)*delta1     ## We generate a random new position for the randomly chosen particle
        y = particulas1[d][1] + (np.random.rand()-0.5)*delta1


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas1[d]=[x,y]            ## If it's a possible step, we accept the new configuration
            

      z = particulas1-particulas_i  ## To calculate MSD
       
      s = []              ## We reset s, for the next iteration
       
      for m in range(len(z)):    ## We discount the PBC from the movement of the particles, to properly estimate the MSD
           
           if(z[m][0]>(L-sigma/2)/2):
               
               z[m][0] = z[m][0] - (L-sigma/2)
           
           if(z[m][1]>(L-sigma/2)/2):
               
               z[m][1] = z[m][1] - (L-sigma/2)
           
           if(z[m][0]<-(L-sigma/2)/2):
               
               z[m][0] = z[m][0] + L-sigma/2
           
           if(z[m][1]<-(L-sigma/2)/2):
               
               z[m][1] = z[m][1] + L-sigma/2
               
      for k in range(N): ## We calculate the MSD in each iteration
           
       
         MSD1[i] = MSD1[i] + 1/N*(z[k][0]**2+z[k][1]**2)

            

    




## Now the process is repeated for all of the perturbation amplitudes.





delta2 = 0.003           ## The amplitude of the random motion

particulas2 = copy.deepcopy(particulas_i)  ## We save the initial positions


MSD2 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    
      while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*N)           ## We choose a random particle
    
        trial = copy.deepcopy(particulas2)
    
        x = particulas2[d][0] + (np.random.rand()-0.5)*delta2     ## We generate a random new position for the random particle
        y = particulas2[d][1] + (np.random.rand()-0.5)*delta2


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas2[d]=[x,y]
            

      z = particulas2-particulas_i
       
      s = []              ## We reset s
       
      for m in range(len(z)):    ## We discount the PBC from the movement of the particles
           
           if(z[m][0]>(L-sigma/2)/2):
               
               z[m][0] = z[m][0] - (L-sigma/2)
           
           if(z[m][1]>(L-sigma/2)/2):
               
               z[m][1] = z[m][1] - (L-sigma/2)
           
           if(z[m][0]<-(L-sigma/2)/2):
               
               z[m][0] = z[m][0] + L-sigma/2
           
           if(z[m][1]<-(L-sigma/2)/2):
               
               z[m][1] = z[m][1] + L-sigma/2
               
      for k in range(N): ##MSD
           
       
         MSD2[i] = MSD2[i] + 1/N*(z[k][0]**2+z[k][1]**2)









delta3 = 0.01           ## The amplitude of the random motion

particulas3 = copy.deepcopy(particulas_i)  ## We save the initial positions



MSD3 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop

    
    while len(s)<1000:  ## Now the Monte Carlo Simulation begins
    
        d = int(np.random.rand()*N)           ## We choose a random particle
    
        trial = copy.deepcopy(particulas3)
    
        x = particulas3[d][0] + (np.random.rand()-0.5)*delta3     ## We generate a random new position for the random particle
        y = particulas3[d][1] + (np.random.rand()-0.5)*delta3


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas3[d]=[x,y]
            
            
    
    z = particulas3-particulas_i
    
    s = []              
    
    for m in range(len(z)):    ## We discount the PBC from the movement of the particles
        
        if(z[m][0]>(L-sigma/2)/2):
            
            z[m][0] = z[m][0] - (L-sigma/2)
        
        if(z[m][1]>(L-sigma/2)/2):
            
            z[m][1] = z[m][1] - (L-sigma/2)
        
        if(z[m][0]<-(L-sigma/2)/2):
            
            z[m][0] = z[m][0] + L-sigma/2
        
        if(z[m][1]<-(L-sigma/2)/2):
            
            z[m][1] = z[m][1] + L-sigma/2
            
            
    for k in range(N):          
    
            MSD3[i] = MSD3[i] + 1/N*(z[k][0]**2+z[k][1]**2)
            
            
            
            
            
            
    
            
delta4 = 0.03           ## The amplitude of the random motion

particulas4 = copy.deepcopy(particulas_i)  ## We save the initial positions


MSD4 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    

    while len(s)<1000:  ## Now the Monte Carlo Simulation begins
    
        d = int(np.random.rand()*N)           ## We choose a random particle
    
        trial = copy.deepcopy(particulas4)
    
        x = particulas4[d][0] + (np.random.rand()-0.5)*delta4     ## We generate a random new position for the random particle
        y = particulas4[d][1] + (np.random.rand()-0.5)*delta4


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas4[d]=[x,y]


    z = particulas4-particulas_i
    
    s = []              
    
    for m in range(len(z)):    ## We discount the PBC from the movement of the particles
        
        if(z[m][0]>(L-sigma/2)/2):
            
            z[m][0] = z[m][0] - (L-sigma/2)
        
        if(z[m][1]>(L-sigma/2)/2):
            
            z[m][1] = z[m][1] -(L-sigma/2)
        
        if(z[m][0]<-(L-sigma/2)/2):
            
            z[m][0] = z[m][0] + L-sigma/2
        
        if(z[m][1]<-(L-sigma/2)/2):
            
            z[m][1] = z[m][1] + L-sigma/2
            
            
    for k in range(N):          
        

        MSD4[i] = MSD4[i] + 1/N*(z[k][0]**2+z[k][1]**2)
    
    







            

delta5 = 0.1           ## The amplitude of the random motion

particulas5 = copy.deepcopy(particulas_i)  ## We save the initial positions


MSD5 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    
    
    while len(s)<1000:  ## Now the Monte Carlo Simulation begins
    
        d = int(np.random.rand()*N)           ## We choose a random particle
    
        trial = copy.deepcopy(particulas5)
    
        x = particulas5[d][0] + (np.random.rand()-0.5)*delta5     ## We generate a random new position for the random particle
        y = particulas5[d][1] + (np.random.rand()-0.5)*delta5


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas5[d]=[x,y]
            
    z = particulas5-particulas_i
    
    s = []              
    
    for m in range(len(z)):    ## We discount the PBC from the movement of the particles
        
        if(z[m][0]>(L-sigma/2)/2):
            
            z[m][0] = z[m][0] - (L-sigma/2)
        
        if(z[m][1]>(L-sigma/2)/2):
            
            z[m][1] = z[m][1] - (L-sigma/2)
        
        if(z[m][0]<-(L-sigma/2)/2):
            
            z[m][0] = z[m][0] + L-sigma/2
        
        if(z[m][1]<-(L-sigma/2)/2):
            
            z[m][1] = z[m][1] + L-sigma/2
            
            
    for k in range(N):          
        
        MSD5[i] = MSD5[i] + 1/N*(z[k][0]**2+z[k][1]**2)            
            








delta6 = 0.3           ## The amplitude of the random motion

particulas6 = copy.deepcopy(particulas_i)  ## We save the initial positions


MSD6 = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    

    
    
    while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*N)           ## We choose a random particle
    
        trial = copy.deepcopy(particulas6)
    
        x = particulas6[d][0] + (np.random.rand()-0.5)*delta6     ## We generate a random new position for the random particle
        y = particulas6[d][1] + (np.random.rand()-0.5)*delta6


        if (x>(L/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(L-sigma/2)
        
        if (y>(L/2-sigma/2)):
            
            y = y-(L-sigma/2)
            
        if (x<-(L/2-sigma/2)):
            
            x = x+(L-sigma/2)
            
        if (y<-(L/2-sigma/2)):
            
            y = y+(L-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            particulas6[d]=[x,y]
            
    z = particulas6-particulas_i
    
    s = []              ## We set an empty array, to cut the "while" loop once a Monte Carlo step has been completed
    
    for m in range(len(z)):    ## We discount the PBC from the movement of the particles
        
        if(z[m][0]>(L-sigma/2)/2):
            
            z[m][0] = z[m][0] - (L-sigma/2)
        
        if(z[m][1]>(L-sigma/2)/2):
            
            z[m][1] = z[m][1] - (L-sigma/2)
        
        if(z[m][0]<-(L-sigma/2)/2):
            
            z[m][0] = z[m][0] + L-sigma/2
        
        if(z[m][1]<-(L-sigma/2)/2):
            
            z[m][1] = z[m][1] + L-sigma/2
            
            
    for k in range(N):          ## for/if loop calculates the MSD, t=0, the MSD = 0, so it doesn't give error in the calculations
        
        MSD6[i] = MSD6[i] + (1/N)*((z[k][0])**2+(z[k][1])**2)



## We plot all the MSD in the same graphic, with a logarithmic scale

plt.figure(figsize=(8,6))
plt.xscale("log")
plt.yscale("log")
plt.title('Mean Square Deviation')
plt.plot(t, MSD1, label='δ=0.001')
plt.plot(t, MSD2, label='δ=0.003')
plt.plot(t, MSD3, label='δ=0.01')
plt.plot(t, MSD4, label='δ=0.03')
plt.plot(t, MSD5, label='δ=0.1')
plt.plot(t, MSD6, label='δ=0.3')


plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.legend(loc='best')
plt.show()

## The diffusivity

D = [MSD1[T-1]/(4*T), MSD2[T-1]/(4*T), MSD3[T-1]/(4*T), MSD4[T-1]/(4*T), MSD5[T-1]/(4*T), MSD6[T-1]/(4*T)]
delta = [delta1, delta2, delta3, delta4, delta5, delta6]



plt.figure(figsize=(8,6))           ## We plot the diffusivity
plt.plot(delta, D, '--ro')
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.title('Diffusivity vs noise amplitude')


plt.figure(figsize=(8,6))           ## We plot the square root of the diffusivity, to check if it's a perfect quadratic relation
plt.plot(delta, np.sqrt(D), '--ro')
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.title('Square root of the diffusivity vs noise amplitude')


#particulas = np.array(particulas)
#plt.plot(particulas[:, 0], particulas[:, 1], 'ok')
#plt.plot(particulas_i[:,0], particulas_i[:,1], 'or')
#%%  2)



triangular_lattice = []  ## We set the triangular lattice initially as an empty array

n = 240         ## Number of particles

phi = 0.05      ## Space occupied by the particles/Space occupied by the box

l = np.sqrt(np.pi*n/(4*phi))        ## Side of the box



for i in range(8):     ## We set the triangular lattice so the particles are as close as possible, but not overlapping

    for k in range(15):
        
        if (i==12):
            
            x1 = -7.5 + k
            y1 = -5.625 + 2*i*np.sqrt(0.75)  ## Square root of 0.75 cause of Pitagoras, since the diameter of each one is equal to 1, so all particles are touching, but not overlapping
            

            triangular_lattice.append([x1,y1])
        else:
            
            x1 = -7.5 + k
            y1 = -5.625 + 2*i*np.sqrt(0.75)
            
            x2 = -7 + k 
            y2 = -5.625+(2*i+1)*np.sqrt(0.75)
            
    
    
            triangular_lattice.append([x1,y1])
            triangular_lattice.append([x2,y2])
        


triangular_lattice=np.array(triangular_lattice)

positions = copy.deepcopy(triangular_lattice)  ## These are the positions that will change in each step, to save the initial configuration

 
plt.figure(figsize=(8,6))  ## We plot the initial configuration to check everything's in order
plt.plot(triangular_lattice[:,0], triangular_lattice[:,1], 'ro', markersize=10)
plt.xlim([-l/2,l/2])
plt.ylim([-l/2,l/2])


dt = 1              ## Time passed on each Monte Carlo step
T = 3000*dt              ## Total number of Monter Carlo steps in my simulation

t = np.linspace(0, T, T)  ## Time array
s = []
delta = 0.3             ## Amplitude of perturbation

MSD = np.zeros(len(t))  ## We set the MSD, just like before


for i in range(len(t)):         ## The simulation begins


        

    while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*n)           ## We generate a random number to pick a random particle
    
        trial = copy.deepcopy(positions)   ## We set a trial, to check if the new position verifies all the conditions, while leaving the original configuration untouched
        
        x = positions[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the randomly chosen particle
        y = positions[d][1] + (np.random.rand()-0.5)*delta


        if (x>(l/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(l/2-sigma/2)
        
        if (y>(l/2-sigma/2)):
            
            y = y-(l/2-sigma/2)
            
        if (x<-(l/2-sigma/2)):
            
            x = x+(l/2-sigma/2)
            
        if (y<-(l/2-sigma/2)):
            
            y = y+(l/2-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            positions[d]=[x,y]   ## We accept the new configuration if it doesn't cause any overlapping
            
            
            
    z = positions-triangular_lattice  ## The difference between initial and final positions, to calculate MSD
    
    s = []              ## We reset s, for the next loop
    

    for m in range(len(z)):
        
        if(z[m][0]>(l-sigma/2)/2):    ## We discount the jump from the MSD in case any particle jumps to the other side of the box because of periodic boundary conditions
            
            z[m][0] = z[m][0] - (l-sigma/2)
        
        if(z[m][1]>(l-sigma/2)/2):
            
            z[m][1] = z[m][1] - (l-sigma/2)
        
        if(z[m][0]<-(l-sigma/2)/2):
            
            z[m][0] = z[m][0] + l-sigma/2
        
        if(z[m][1]<-(l-sigma/2)/2):
            
            z[m][1] = z[m][1] + l-sigma/2
                        
    for k in range(n):  ## MSD

        
        MSD[i] = MSD[i] + 1/n*(z[k][0]**2+z[k][1]**2)
        
    if (i==100 or i==200 or i==500 or i==600 or i==800 or i==1000 or i==1200 or i==1500 or i==1700 or i==1800 or i==2999):
        
        ## We plot the system periodically to see how it evolves in time
    
        plt.figure(figsize=[8,6])
        plt.title("ϕ = 0.05, n = %i" %i)
        plt.plot(positions[:,0], positions[:,1], 'or')
        plt.xlim([-l/2,l/2])
        plt.ylim([-l/2,l/2])
        plt.savefig("phi005_" + str(i) + ".png")
        plt.show()
    
        
plt.figure(figsize=[8,6])       ## We plot the MSD to see when it reaches equilibrium
plt.plot(t, MSD, 'r')
plt.xscale("log")
plt.yscale("log") 
plt.title("MSD ϕ = 0.05")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("t")
plt.ylabel("MSD")
plt.savefig("MSD_phi005.png")
plt.show()              
        


## Now we change the side of the box and repeat the process


dt = 1              ## Time passed on each Monte Carlo step
T = 2000*dt              ## Total number of Monter Carlo steps in my simulation

t = np.linspace(0, T, T)   ## Time array
s = []
delta = 0.3             ## Amplitude of perturbation

MSD2_2 = np.zeros(len(t))

positions2 = copy.deepcopy(triangular_lattice)

phi2 = 0.2  ## Relation between the volume occupied by the particles and the total volume of the box

l2 = np.sqrt(np.pi*n/(4*phi2))  ## Length of the side of the new box

for i in range(len(t)):


        

    while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*n)           ## We choose a random particle
    
        trial = copy.deepcopy(positions2)
    
        x = positions2[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = positions2[d][1] + (np.random.rand()-0.5)*delta


        if (x>(l2/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(l2/2-sigma/2)
        
        if (y>(l2/2-sigma/2)):
            
            y = y-(l2/2-sigma/2)
            
        if (x<-(l2/2-sigma/2)):
            
            x = x+(l2/2-sigma/2)
            
        if (y<-(l2/2-sigma/2)):
            
            y = y+(l2/2-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            positions2[d]=[x,y]
            
            
            
    z = positions2-triangular_lattice
    
    s = []  ## We reset s
           
    
    for m in range(len(z)):   ## We discount the jump from the MSD in case any particle jumps to the other side of the box because of periodic boundary conditions
        
        if(z[m][0]>(l2-sigma/2)/2):
            
            z[m][0] = z[m][0] - (l2-sigma/2)
        
        if(z[m][1]>(l2-sigma/2)/2):
            
            z[m][1] = z[m][1] - (l2-sigma/2)
        
        if(z[m][0]<-(l2-sigma/2)/2):
            
            z[m][0] = z[m][0] + (l2-sigma/2)
        
        if(z[m][1]<-(l2-sigma/2)/2):
            
            z[m][1] = z[m][1] + (l2-sigma/2)
                

                        
    for k in range(n):         

        
        MSD2_2[i] = MSD2_2[i] + 1/n*(z[k][0]**2+z[k][1]**2)
        
    if (i==100 or i==200 or i==500 or i==650 or i==800 or i==1000 or i==1200 or i==1500 or i==1700 or i==1800 or i==1999):
        
        ## We periodically plot the system, to check its evolution
        
        plt.figure(figsize=[8,6])
        plt.title("ϕ = 0.2 n = %i" %i)
        plt.plot(positions2[:,0], positions2[:,1], 'or')
        plt.xlim([-l2/2,l2/2])
        plt.ylim([-l2/2,l2/2])
        plt.savefig("phi02_" + str(i) + ".png")
        plt.show()


        
plt.figure(figsize=[8,6])
plt.title("MSD ϕ = 0.2")    
plt.plot(t, MSD2_2, 'r')
plt.xscale("log")
plt.yscale("log") 
plt.xlabel("t")
plt.ylabel("MSD")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.savefig("MSD_phi02.png")
plt.show()


## Now, for the third configuration


MSD3_2 = np.zeros(len(t))

positions3 = copy.deepcopy(triangular_lattice)

phi3 = 0.5

l3 = np.sqrt(np.pi*n/(4*phi3))  

for i in range(len(t)):


        

    while len(s)<1000:  ## Now the Monte Carlo Simulation begins
    
        d = int(np.random.rand()*n)           ## We choose a random particle
    
        trial = copy.deepcopy(positions3)
    
        x = positions3[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = positions3[d][1] + (np.random.rand()-0.5)*delta


        if (x>(l3/2-sigma/2)):           ## Periodidc boundary conditions
            
            x = x-(l3/2-sigma/2)
        
        if (y>(l3/2-sigma/2)):
            
            y = y-(l3/2-sigma/2)
            
        if (x<-(l3/2-sigma/2)):
            
            x = x+(l3/2-sigma/2)
            
        if (y<-(l3/2-sigma/2)):
            
            y = y+(l3/2-sigma/2)
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.

        if check(x, y, trial, diametro=sigma):  ## We check if the new position is acceptable with our conditions
            
            positions3[d]=[x,y]
            
            
            
    z = positions3-triangular_lattice
    
    s = []  ## We reset s
            
    if(z[m][0]>(l3-sigma/2)/2):    ## Discount of PBC from z
            
        z[m][0] = z[m][0] - (l3-sigma/2)
        
    if(z[m][1]>(l3-sigma/2)/2):
            
        z[m][1] = z[m][1] - (l3-sigma/2)
        
    if(z[m][0]<-(l3-sigma/2)/2):
            
        z[m][0] = z[m][0] + (l3-sigma/2)
        
    if(z[m][1]<-(l3-sigma/2)/2):
            
        z[m][1] = z[m][1] + (l3-sigma/2)
                

                        
    for k in range(n):          ## MSD

        
        MSD3_2[i] = MSD3_2[i] + 1/n*(z[k][0]**2+z[k][1]**2)
        
    if (i==100 or i==200 or i==500 or i==650 or i==800 or i==1000 or i==1200 or i==1500 or i==1700 or i==1800 or i==1999):
        
        
        plt.figure(figsize=[8,6])
        plt.title("ϕ = 0.5 n = %i" %i)
        plt.plot(positions3[:,0], positions3[:,1], 'or')
        plt.xlim([-l3/2,l3/2])
        plt.ylim([-l3/2,l3/2])
        plt.savefig("phi05_" + str(i) + "_.png")
        plt.show()
        
plt.figure(figsize=[8,6])
plt.title("MSD ϕ = 0.5")
plt.xscale("log")
plt.yscale("log")        
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.plot(t, MSD3_2, 'r')
plt.xlabel("t")
plt.ylabel("MSD")
plt.savefig("MSD_phi05.png")
plt.show()

#%% 5)


Lx = 15

Ly = 10*Lx      ## We set the box's size

g1 = 0.01       ## Gravity 

Temp = 1        ## Temperature
m = 1           ## Mass of the particles
kb = 1            ## Boltzmann's constant
N = 1000            ## Number of particles
sigma = 1


MV = np.zeros(5)   ## We set an empty array to save the mean value of the y component at the end of the Monte Carlo Simulations

## We use the previously defined check function, which verifies that the random particle 
#  generated does not overlap with the rest of the particles.



def check(x, y, particulas, diametro=1):   
    cond = all(
        [(np.sqrt((x - p[0]) ** 2 + (y - p[1]) ** 2) >= sigma) for p in particulas]
    )
    return cond

np.random.seed()

r_i = []                 ## Array of particles


## We create a while loop that creates particles randomly distributed within the box
#  and we include the particles who do not overlap to our system and discard
#  those that overlap.

while len(r_i) < 1000:  ## We set the initial configuration

    x_i = np.random.rand() * (Lx - sigma) + sigma/2
    y_i = np.random.rand() * (Ly - sigma) + sigma/2

    if check(x_i, y_i, r_i, diametro=sigma):
        r_i.append([x_i, y_i])

## We plot the particles to check that they, indeed, do not overlap


r_i = np.array(r_i)






plt.figure(figsize=(8,16))
plt.plot(r_i[:, 0], r_i[:, 1], "or")
plt.xlim([0,Lx])
plt.ylim([0,Ly])
plt.title("g = 0")
plt.xlabel("x")
plt.ylabel("y")   
plt.show()
  


delta = 0.3            ## The amplitude of the random motion, it will always be 0.3

s = []      ## We set an empty array to cut the while loop, like before

g0 = 0  ## Value of the gravity

r0 = copy.deepcopy(r_i)  ## The positions that will evolve with each iteration

dt = 1               ## Time jump, the unit of time is one Monte Carlo step

T = 1000*dt


t = np.linspace(0, T, int(T/dt))        ## Our time

MSD0 = np.zeros(len(t))   ## We set the MSD

for i in range(len(t)):     ## We begin the simulation loop
    

    
    while len(s)<1000:  ## Now a Monte Carlo Step begins
    
        d = int(np.random.rand()*1000)           ## We generate a random number to choose a particle randomly
    
        trial = copy.deepcopy(r0)
        
        y_i = r0[d][1]                   ## We save the initial y component of the random particle selected, to calculate the probability of transition
        
        x = r0[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = r0[d][1] + (np.random.rand()-0.5)*delta


        if (x > (Lx-sigma/2)):           ## Hard wall interaction condition, if the particle tries to go further than the wall, it will "collide" with it
            
            x = Lx-sigma/2
        
        if (y > (Ly-sigma/2)):
            
            y = Ly-sigma/2
            
        if (x < sigma/2):
            
            x = sigma/2
            
        if (y < sigma/2):
            
            y = sigma/2
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.
        
        h = np.random.rand()     ## We generate a random number, to see if we accept the new configuration, based on the probability of transition
        
        if (check(x, y, trial, diametro=sigma) and (h<=(np.exp(-1/(kb*Temp)*m*g0*(y-y_i))))):  ## We check if the new position is acceptable with our conditions, and if the generated random number is smaller than the probability of transiton
            
            r0[d]=[x,y]   ## If the new position is acceptable and the random nubmer generated is smaller than the probability of transition, the new configuration is accepted
            
            
            
    z = r0-r_i
    
    s = []              
    
            
            
    for k in range(N):          ## for/if loop calculates the MSD, t=0, the MSD = 0, so it doesn't give error in the calculations

        
        MSD0[i] = MSD0[i] + 1/N*(z[k][0]**2+z[k][1]**2)
    

MV[0] = np.sum(r0[:,1])/N  ## The mean value of y


plt.figure(figsize=[8,16])      ##We plot the final configuration
plt.plot(r0[:,0], r0[:,1], 'or', label = "Final configuration")
#plt.plot(r_i[:,0], r_i[:,1], 'ok', label = "Initial configuration")    ##Uncomment if you want a comparative between the final and the initial positions
#plt.legend(loc="upper right")
plt.xlim([0,Lx])
plt.ylim([0,Ly])
plt.title("g = 0")
plt.xlabel("x")
plt.ylabel("y")     

plt.savefig("g0.png")
plt.show()



plt.figure(figsize=[8,6])  ## We plot the MSD, to verify equilibrium conditions
plt.plot(t, MSD0, 'r')
plt.title("MSD g = 0.01")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("t")
plt.ylabel("MSD")
plt.xscale("log")
plt.yscale("log") 
plt.savefig("MSD_g_0.png")

plt.show()


## We repeat the process for the different values of g

dt = 1               ## Time jump, the unit of time is one Monte Carlo step

T = 1000*dt

s = []

r = copy.deepcopy(r_i)

t = np.linspace(0, T, int(T/dt))        ## Our time

MSD_ = np.zeros(len(t))          ## We set the Mean Square Deviation

for i in range(len(t)):     ## We begin the simulation loop
    

    
    while len(s)<1000:  ## Now the Monte Carlo Step begins
    
        d = int(np.random.rand()*1000)           ## We choose a random particle
    
        trial = copy.deepcopy(r)
        
        y_i = r[d][1]
        
        x = r[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = r[d][1] + (np.random.rand()-0.5)*delta


        if (x > (Lx-sigma/2)):           ## Hard wall interaction condition
            
            x = Lx-sigma/2
        
        if (y > (Ly-sigma/2)):
            
            y = Ly-sigma/2
            
        if (x < sigma/2):
            
            x = sigma/2
            
        if (y < sigma/2):
            
            y = sigma/2
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.
        
        h = np.random.rand()
        
        if (check(x, y, trial, diametro=sigma) and (h<=(np.exp(-1/(kb*Temp)*m*g1*(y-y_i))))):  ## We check if the new position is acceptable with our conditions
            
            r[d]=[x,y]
            
            
            
    z = r-r_i
    
    s = []              
    
            
            
    for k in range(N):         

        
        MSD_[i] = MSD_[i] + 1/N*(z[k][0]**2+z[k][1]**2)
    

   
plt.figure(figsize=[8,16])
plt.plot(r[:,0], r[:,1], 'or', label = "Final configuration")
#plt.plot(r_i[:,0], r_i[:,1], 'ok', label = "Initial configuration")
#plt.legend(loc="upper right"))
plt.xlim([0,Lx])
plt.ylim([0,Ly])
plt.title("g = 0.01")
plt.xlabel("x")
plt.ylabel("y")     

plt.savefig("g0_01.png")
plt.show()



plt.figure(figsize=[8,6])
plt.plot(t, MSD_, 'r')
plt.title("MSD g = 0.01")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("t")
plt.ylabel("MSD")
plt.xscale("log")
plt.yscale("log") 
plt.savefig("MSD_g_0_01.png")

plt.show()


MV[1] = np.sum(r[:,1]/N)
        
    
        


## Now for g = 0.1

g2 = 0.1

r2 = copy.deepcopy(r_i)
dt = 1              
T = 1000*dt
t = np.linspace(0, T, int(T/dt))
MSD3_2 = np.zeros(len(t))

for i in range(len(t)):     ## We begin the simulation loop
    

    
    while len(s)<1000:  ## Now a Monte Carlo Step begins
    
        d = int(np.random.rand()*1000)           ## We choose a random particle
    
        trial = copy.deepcopy(r2)
        
        y_i = r2[d][1]
        
        x = r2[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = r2[d][1] + (np.random.rand()-0.5)*delta


        if (x > (Lx-sigma/2)):           ## Hard wall interaction condition
            
            x = Lx-sigma/2
        
        if (y > (Ly-sigma/2)):
            
            y = Ly-sigma/2
            
        if (x < sigma/2):
            
            x = sigma/2
            
        if (y < sigma/2):
            
            y = sigma/2
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.
        
        h = np.random.rand()
        
        if (check(x, y, trial, diametro=sigma) and (h<=(np.exp(-1/(kb*Temp)*m*g2*(y-y_i))))):  ## We check if the new position is acceptable with our conditions
            
            r2[d]=[x,y]
            
            
            
    z = r2-r_i
    
    s = []              
    
            
            
    for k in range(N):         

        
        MSD3_2[i] = MSD3_2[i] + 1/N*(z[k][0]**2+z[k][1]**2)
    

        
        
        
plt.figure(figsize=[8,16])
plt.plot(r2[:,0], r2[:,1], 'or', label = "Final configuration")
#plt.plot(r_i[:,0], r_i[:,1], 'ok', label = "Initial configuration")
#plt.legend(loc="upper right"))
plt.title("g = 0.1")
plt.xlim([0,Lx])
plt.ylim([0,Ly])
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("g_0_1.png")
plt.show()



plt.figure(figsize=[8,6])
plt.plot(t, MSD3_2, 'r')
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.title("MSD_g_01")
plt.xlabel("t")
plt.xscale("log")
plt.yscale("log") 
plt.ylabel("MSD")
plt.savefig("MSD_g_0_1.png")
plt.show()

print(np.sum(r_i[:,1]/N))
MV[2] = np.sum(r2[:,1]/N)
    
        
    
    
    
## For g = 1
    
    
g3 = 1

r3 = copy.deepcopy(r_i)
dt = 1              
T = 2000*dt
t = np.linspace(0, T, int(T/dt))
MSD3_3 = np.zeros(len(t))


for i in range(len(t)):     ## We begin the simulation loop
    

    
    while len(s)<1000:  ## Now the Monte Carlo Simulation begins
    
        d = int(np.random.rand()*1000)           ## We choose a random particle
    
        trial = copy.deepcopy(r3)
        
        y_i = r3[d][1]
        
        x = r3[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = r3[d][1] + (np.random.rand()-0.5)*delta


        if (x > (Lx-sigma/2)):           ## Hard wall interaction condition
            
            x = Lx-sigma/2
        
        if (y > (Ly-sigma/2)):
            
            y = Ly-sigma/2
            
        if (x < sigma/2):
            
            x = sigma/2
            
        if (y < sigma/2):
            
            y = sigma/2
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.
        
        h = np.random.rand()
        
        if (check(x, y, trial, diametro=sigma) and (h<=(np.exp(-1/(kb*Temp)*m*g3*(y-y_i))))):  ## We check if the new position is acceptable with our conditions
            
            r3[d]=[x,y]
            
            
            
    z = r3-r_i
    
    s = []             
    
            
            
    for k in range(N):         

        
        MSD3_3[i] = MSD3_3[i] + 1/N*(z[k][0]**2+z[k][1]**2)
        
        
    if (i==100 or i==200 or i==500 or i==650 or i==800 or i==1000 or i==1200 or i==1499):




       
        plt.figure(figsize=[8,16])
        plt.title("g = 1")
        plt.plot(r3[:,0], r3[:,1], 'or', label = "Final configuration")
        #plt.plot(r_i[:,0], r_i[:,1], 'ok', label = "Initial configuration")
        #plt.legend(loc="upper right"))
        plt.xlim([0,Lx])
        plt.ylim([0,Ly])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig("g_1_"+ str(i) +"_.png")
        plt.show()
        
        
MV[3] = np.sum(r3[:,1]/N)


plt.figure(figsize=[8,6])
plt.plot(t, MSD3_3, 'r')
plt.title("MSD g = 1")
plt.xscale("log")
plt.yscale("log") 
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("t")
plt.ylabel("MSD")
plt.savefig('MSD_g_1.png')
plt.show()        
       
        
       
## And g = 10

g4 = 10

r4 = copy.deepcopy(r_i)
dt = 1              
T = 3000*dt

t = np.linspace(0, T, int(T/dt))
MSD3_4 = np.zeros(len(t))

for i in range(len(t)):     ## We begin the simulation loop
    

    
    while len(s)<1000:  ## Now a Monte Carlo Step begins
    
        d = int(np.random.rand()*1000)           ## We choose a random particle
    
        trial = copy.deepcopy(r4)
        
        y_i = r4[d][1]
        
        x = r4[d][0] + (np.random.rand()-0.5)*delta     ## We generate a random new position for the random particle
        y = r4[d][1] + (np.random.rand()-0.5)*delta


        if (x > (Lx-sigma/2)):           ## Hard wall interaction condition
            
            x = Lx-sigma/2
        
        if (y > (Ly-sigma/2)):
            
            y = Ly-sigma/2
            
        if (x < sigma/2):
            
            x = sigma/2
            
        if (y < sigma/2):
            
            y = sigma/2
            
            
        s.append(1)
        
        trial[d] = [999999999,999999999]        ## We set the position of the moved particle far away, since otherwise it would be overlaping with the new possible configuration and the new position would mistakenly never be acceptable.
        
        h = np.random.rand()
        
        if (check(x, y, trial, diametro=sigma) and (h<=(np.exp(-1/(kb*Temp)*m*g4*(y-y_i))))):  ## We check if the new position is acceptable with our conditions
            
            r4[d]=[x,y]
            
            
            
    z = r4-r_i
    
    s = []              
    
            
            
    for k in range(N):          

        
        MSD3_4[i] = MSD3_4[i] + 1/N*(z[k][0]**2+z[k][1]**2)
    
    
    if (i==100 or i==200 or i==500 or i==650 or i==800 or i==1000 or i==1200 or i==1500 or i==1700 or i==1900 or i==1999 or i==2300 or i==2500 or i==2700 or i==2900 or i==2999):
        
        
        plt.figure(figsize=[8,16])
        plt.title("g = 10")
        plt.plot(r4[:,0], r4[:,1], 'or', label = "Final configuration")
        #plt.plot(r_i[:,0], r_i[:,1], 'ok', label = "Initial configuration")
        #plt.legend(loc="upper right"))
        plt.xlim([0,Lx])
        plt.ylim([0/2,Ly])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig("g10_" + str(i) + "_.png")
        plt.show()


MV[4] = np.sum(r4[:,1])/N




plt.figure(figsize=[8,16])
plt.plot(r4[:,0], r4[:,1], 'or', label="final configuration")
plt.plot(r_i[:,0], r_i[:,1], 'ok', label="initial configuration")
plt.legend(loc="best")
plt.xlim([0,Lx])
plt.ylim([0,Ly])
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("g_10.png")
plt.show()



plt.figure(figsize=[8,6])
plt.plot(t, MSD3_4, 'r')
plt.xscale("log")
plt.yscale("log") 
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.xlabel("t")
plt.ylabel("MSD")
plt.savefig("MSD_g_10.png")
plt.show()
    
## Now we write the mean values of y in a .txt
        
file = open("C:/Users/Usuario/Desktop/Project_2/Apartado_5/MV.txt", "w+")  

for i in range(len(MV)):
    content = str(MV[i])
    file.write("g = %i" %i)
    file.write('\n')
    file.write(content)
    file.write('\n')
    
file.close()

    
    



