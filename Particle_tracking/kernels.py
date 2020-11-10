import math


def SampleTemp(particle, fieldset, time):
    particle.temp = fieldset.T[time, particle.depth, particle.lat, particle.lon]
    particle.salin = fieldset.S[time, particle.depth, particle.lat, particle.lon]

def Sink(particle, fieldset, time):
    if particle.depth > fieldset.dwellingdepth:
        particle.depth = particle.depth + fieldset.sinkspeed * particle.dt
    else:
        particle.depth = fieldset.dwellingdepth

def initial(particle, fieldset, time):
    if(not (particle.age>-1)): # if particle.age==nan, then put it on the bathymetry
        particle.age = 0
        particle.depth = fieldset.B[time, fieldset.surface, particle.lat, particle.lon]
        if(particle.depth<fieldset.dwellingdepth): # Directly delete the particle if it is below the dwelling depth
            particle.delete()

def Age(particle, fieldset, time):
    if particle.depth <= fieldset.dwellingdepth:
        particle.age = particle.age + math.fabs(particle.dt)
    if particle.age > fieldset.maxage:
        particle.delete()

def DeleteParticle(particle, fieldset, time):
    particle.delete()

#%% Kernels that use mixed layer depth as dwelling depth
def initialML(particle, fieldset, time):
    if(not (particle.age>-1)): # if particle.age==nan, then put it on the bathymetry
        particle.age = 0
        particle.depth = fieldset.B[time, fieldset.surface, particle.lat, particle.lon]

def SinkML(particle, fieldset, time):
    if(particle.MLr==0): # sinking if particle did not reach dwelling depth yet
        particle.depth = particle.depth + fieldset.sinkspeed * particle.dt
        if(particle.depth < fieldset.ML[time, particle.depth, particle.lat, particle.lon]):
            particle.MLr = 1
    else: # put particle on mixed layer depth otherwise
        particle.depth = fieldset.ML[time, particle.depth, particle.lat, particle.lon]

def AgeML(particle, fieldset, time):
    if particle.MLr==1:
        particle.age = particle.age + math.fabs(particle.dt)
    if particle.age > fieldset.maxage:
        particle.delete()
