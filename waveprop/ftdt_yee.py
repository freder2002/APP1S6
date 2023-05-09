# GPL3, Copyright (c) Max Hofheinz, UdeS, 2021

import numpy, fiddle
from subprocess import Popen, PIPE
import mmap
import time


def curl_E(E):
    curl_E = numpy.zeros(E.shape)
    curl_E[:, :-1, :, 0] += E[:, 1:, :, 2] - E[:, :-1, :, 2]
    curl_E[:, :, :-1, 0] -= E[:, :, 1:, 1] - E[:, :, :-1, 1]

    curl_E[:, :, :-1, 1] += E[:, :, 1:, 0] - E[:, :, :-1, 0]
    curl_E[:-1, :, :, 1] -= E[1:, :, :, 2] - E[:-1, :, :, 2]

    curl_E[:-1, :, :, 2] += E[1:, :, :, 1] - E[:-1, :, :, 1]
    curl_E[:, :-1, :, 2] -= E[:, 1:, :, 0] - E[:, :-1, :, 0]
    return curl_E

def curl_H(H):
    curl_H = numpy.zeros(H.shape)

    curl_H[:,1:,:,0] += H[:,1:,:,2] - H[:,:-1,:,2]
    curl_H[:,:,1:,0] -= H[:,:,1:,1] - H[:,:,:-1,1]

    curl_H[:,:,1:,1] += H[:,:,1:,0] - H[:,:,:-1,0]
    curl_H[1:,:,:,1] -= H[1:,:,:,2] - H[:-1,:,:,2]

    curl_H[1:,:,:,2] += H[1:,:,:,1] - H[:-1,:,:,1]
    curl_H[:,1:,:,2] -= H[:,1:,:,0] - H[:,:-1,:,0]
    return curl_H

def subp():
    subproc = Popen(["./app", FNAME], stdin=PIPE, stdout=PIPE)
    return subproc

def signal_and_wait(subproc):
    subproc.stdin.write("START\n".encode())
    subproc.stdin.flush()                   # Nécessaire pour vider le tampon de sortie
    res = subproc.stdout.readline()
    #print(res)

FNAME       = "APP-shm"
MATRIX_SIZE = 100

# Lancement de l'exécutable associé
# NOTE: suppose que l'exécutable est dans le même dossier que celui en cours (normalement build/)
subproc = subp()

# Envoi d'une ligne sur l'entrée du sous-processus et attend un retour pour signaler que
# nous sommes prêts à passer à la prochaine étape. 
signal_and_wait(subproc)

shm_f = open(FNAME, "r+b")
shm_mm = mmap.mmap(shm_f.fileno(), 0)

def timestep(E, H, courant_number, source_pos, source_val):
    # methode

    shared_matrix = numpy.ndarray(shape=(2,3,MATRIX_SIZE, MATRIX_SIZE,MATRIX_SIZE), dtype=numpy.float64, buffer=shm_mm)
    # print(shared_matrix)

    # print(source_pos)
    shared_matrix[0,0,0,0,0] = courant_number
    shared_matrix[0,0,0,0,1] = source_pos[0][0]
    shared_matrix[0,0,0,0,2] = source_pos[1][0]
    shared_matrix[0,0,0,0,3] = source_pos[2][0]
    shared_matrix[0,0,0,0,4] = source_pos[3][0]
    shared_matrix[0,0,0,0,5] = source_val  


    signal_and_wait(subproc)
    E = shared_matrix[0,:,:,:,:]
    H = shared_matrix[1,:,:,:,:]

    # print(E)
    # print(H)
    # ancienne methode 
    # E += courant_number * curl_H(H)
    # E[source_pos] += source_val
    # H -= courant_number * curl_E(E)
    return E, H


class WaveEquation:
    def __init__(self, s, courant_number, source):
        s = s + (3,)
        self.E = numpy.zeros(s)
        self.H = numpy.zeros(s)
        self.courant_number = courant_number
        self.source = source
        self.index = 0

    def __call__(self, figure, field_component, slice, slice_index, initial=False):
        if field_component < 3:
            field = self.E
        else:
            field = self.H
            field_component = field_component % 3
        if slice == 0:
            field = field[slice_index, :, :, field_component]
        elif slice == 1:
            field = field[:, slice_index, :, field_component]
        elif slice == 2:
            field = field[:, :, slice_index, field_component]
        source_pos, source_index = source(self.index)
        self.E, self.H = timestep(self.E, self.H, self.courant_number, source_pos, source_index)

        if initial:
            axes = figure.add_subplot(111)
            self.image = axes.imshow(field, vmin=-1e-2, vmax=1e-2)
        else:
            self.image.set_data(field)
        self.index += 1


if __name__ == "__main__":
    n = 100
    r = 0.01
    l = 30

    # Lancement de l'exécutable associé
    # NOTE: suppose que l'exécutable est dans le même dossier que celui en cours (normalement build/)
    subproc = subp()

    # initial wait
    signal_and_wait(subproc)


    def source(index):
        return ([n // 3], [n // 3], [n // 2],[0]), 0.1*numpy.sin(0.1 * index)


    w = WaveEquation((n, n, n), 0.1, source)
    fiddle.fiddle(w, [('field component',{'Ex':0,'Ey':1,'Ez':2, 'Hx':3,'Hy':4,'Hz':5}),('slice',{'XY':2,'YZ':0,'XZ':1}),('slice index',0,n-1,n//2,1)], update_interval=0.01)
