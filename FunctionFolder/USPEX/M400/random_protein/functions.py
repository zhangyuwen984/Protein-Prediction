#!/usr/bin/env python
# encoding: utf-8

'''
Created on Feb 05, 2014

@author: mrakitin
'''

import os, sys
from random import randint


#------------------------------------------------------------------------------#
#                                  Functions:                                  #
#------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------
# Generate angles randomly from range [-180, 180]
def unlimited_random_angles():
    limit_min = -180
    limit_max =  180

    phi = randint(limit_min, limit_max)
    psi = randint(limit_min, limit_max)

    random_angles = {
                     'phi': phi,
                     'psi': psi,
                    }

    return random_angles
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Get random phi and psi pair out of 7 known values for proteins:
def pseudo_random_angles():

    real_angles = [
               {
                'name': 'Right-handed alpha-helix',
                'id'  :   1,
                'phi' : -57,
                'psi' : -47,
               }, {
                'name': '3_10 helix',
                'id'  :   2,
                'phi' : -49,
                'psi' : -26,
               }, {
                'name': 'Parallel beta-sheet',
                'id'  :   3,
                'phi' : -119,
                'psi' :  113,
               }, {
                'name': 'Antiparallel beta-sheet',
                'id'  :   4,
                'phi' : -139,
                'psi' :  135,
               }, {
                'name': 'Left-handed polyproline I',
                'id'  :   5,
                'phi' : -83,
                'psi' : 158,
               }, {
                'name': 'Right-handed polyproline II',
                'id'  :   6,
                'phi' : -78,
                'psi' : 149
               }, {
                'name': 'Pi-helix (rare in protein)',
                'id'  :   7,
                'phi' :  55,
                'psi' : -70
               },
              ]

    phi_psi_id = randint(1, len(real_angles)) - 1

    return real_angles[phi_psi_id]
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Get random angles from the specified list of angles:
def custom_random_angles(custom_angles):

    phi_psi_id = randint(1, len(custom_angles)) - 1

    return custom_angles[phi_psi_id]
#-------------------------------------------------------------------------------

