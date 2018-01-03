# Run FDTD

#Finite Difference Time Domain (FDTD) method for room acoustic simulation
using AcFDTD2
using DSP

#X = 0.01                   # spatial sampling
#env = FDTDEnv(X,IISO())   # create new acoustic env with default values

#where `IISO()` returns the Interpolated Isotropic scheme.
#Alternatively one can choose a samplng frequency 
#instead of a spatial sampling:
#```julia

Fs = 20000.                          # sampling frequency in Hz
env = FDTDEnv(IISO(),Fs; c = 343)   # create new acoustic env with default values

#```
#notice that in the latter line the speed of sound was 
#chosen to be `340` m/s. 
#By default this is set to `343` m/s.
#Set the acoustic impedance `ξ` and room geometry room geometry:
#```julia

ξ = [50.;50.;50.;50.;50.;50.]; # [   ξx1    ;    ξx2   ;    ξy1   ;    ξy2    ;  ξz1 ;   ξz2  ]
                                # [front wall; rear wall; left wall; right wall; floor; ceiling]
geo = CuboidRoom(5, 4, 3, ξ, env)

#```
#The first three parameters 
#indicate the number of spatial samples 
#of the `x`, `y` and `z` directions.  
#Alternatively one can specify 
#the room dimensions in meters:
#```julia
#geo = CuboidRoom(4., 5., 3., ξ, env)
#```
#which are then approximated on the grid.
#Set the number of time steps `Nt`
#Create a band-limited sound source
#with e.g. the `DSP` package:
#```julia


Nt = round(Int64,env.Fs)         # number of time steps (1 sec)
s = zeros(Nt)                    # source signal
s[3] = 1
f2 = geo.env.Fs/2*0.175          # cut-off frequency of source
filt!(s,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(5)),s)

#```
#Define the position of microphone 
#and sound sources:
#```julia

xr = [2 2 2; 2 geo.Ny-1 geo.Nz-1]' # mic positions
xs = [geo.Nx-1 geo.Ny-1 geo.Nz-1]' # sound source position

#must be Array{Int64} with size(xr,1) = 3
#```
#Now type:
#```julia
p = fdtd(s,xs,xr,Nt,geo)
writedlm("FDTD.txt", p)

#```
#to obtain the sound pressure of the microphones.

#For more details on the methods type:
#```julia
#?fdtd
#```

## Credits

#AcFDTD.jl is developed by [Niccolò Antonello](http://homes.esat.kuleuven.be/~nantonel/) at [KU Leuven, ESAT/Stadius](https://www.esat.kuleuven.be/stadius/).











