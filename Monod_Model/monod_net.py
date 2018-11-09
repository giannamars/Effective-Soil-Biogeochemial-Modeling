# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 12:17:17 2018

@author: Gianna
"""

from SloppyCell.ReactionNetworks import *
net = IO.from_SBML_file('Monod.xml', 'monod')

net.set_var_optimizable('Y', False)

pred_net = net.copy()
pred_net.set_var_ic('x_1', net.get_var_ic('x_1')*0.75)

