#!/usr/bin/env python
#
#    Copyright 2011 Gregor Burger
#
#    This file is part of py_epanet.
#
#    py_epanet is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    py_epanet is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with py_epanet.  If not, see <http://www.gnu.org/licenses/>.
#


__author__ = 'Gregor Burger <gregor.burger@uibk.ac.at.at>'

import sys
import input
import math
import collections
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg
import numpy as np

PI = 3.141592654
#PI = math.pi


class UnitsUS:
    def diam(self, dia):
        return dia/12.0

    def length(self, length):
        return length

    def from_gpm(self, gpm):
        return gpm/448.831 # 448.8 = 7.48*60

    def from_hcf(self, hcf):
        return hcf

    def convert_head(self, head):
        return head



class Link:
    def __init__(self, link_data, nodes):
        self.id = link_data.ID
        self.node1 = nodes[link_data.Node1]
        self.node2 = nodes[link_data.Node2]
        self.r = self.resistance() #resistance coefficient from table 3.1 epanet manual
        self.id_headloss = None
        self.flow_correction_factor = None
        self.upper_matrix_offset = None  #direct offset into csr data array
        self.lower_matrix_offset = None  #direct offset into csr data array

    def update(self, Q):
        raise NotImplemented


def power_curve(h0, h1, h2, q1, q2):
    a = h0
    h4 = h0 - h1
    h5 = h0 - h2
    c = math.log(h5/h4)/math.log(q2/q1)

    if c <= 0.0 or c > 20.0:
        raise ValueError

    b = -h4/math.pow(q1, c)

    if b >= 0.0:
        raise ValueError

    return a, b, c

class Pump(Link):
    def __init__(self, pump_data, curves, nodes, units):
        Link.__init__(self, pump_data, nodes)
        self.curve = curves[pump_data.parameters['HEAD']]

        n = len(self.curve.x)

        self.k = 1.0 # somehow strange link configuration K[] in vars.h

        if n == 1:
            self.q1 = self.curve.x[0]#units.from_gpm(self.curve.x[0])
            self.h1 = self.curve.y[0]#units.from_gpm(self.curve.y[0])
            self.h0 = 1.33334*self.h1
            self.q2 = 2.0 * self.q1
            self.h2 = 0.0

            a, b, c = power_curve(self.h0, self.h1, self.h2, self.q1, self.q2)

            self.h0 = -a
            ucf = 448.831 #TODO hack dirty dirty hack use a conversion
            self.R  = -b #* math.pow(ucf, c)
            self.n  = c
            self.q0 = self.q1 #units.from_gpm(self.q1)
            self.q_max  = math.pow((-a/b),(1.0/c))# / ucf
            self.h_max  = -self.h0# / ucf


            #convert units
            self.R *= math.pow(ucf, self.n) # /ucf[HEAD]
            self.q0 /= ucf
            self.q_max /= ucf
            self.h_max /= ucf
        else:
            raise NotImplemented

        self.speed = 1.0
        if 'SPEED' in pump_data.parameters:
            self.speed = float(pump_data.parameters['SPEED'])


    def resistance(self):
        return 1.e8

    def initial_flow(self):
        #print "%30.15f" % self.q0
        #print self.speed * self.q0
        #hydraul.c:324
        return self.speed * self.q0

    def update(self, Q):
        self.id_headloss = self.calc_id_headloss(Q)
        self.flow_correction_factor = self.calc_flow_correction_factor(self.id_headloss, Q)


    def calc_id_headloss(self, Q):
        r = self.R * math.pow(self.k, 2.0-self.n)
        r = self.n*r*math.pow(abs(Q), self.n-1.0)
        return 1.0/max(r, EPAnet.RQtol)

    def calc_flow_correction_factor(self, head_loss, Q):
        h0 = self.k*self.k * self.h0
        n = self.n
        return Q/n + head_loss*h0

    def correct_flow(self):
        #TODO
        """
        Prevent flow in constant HP pumps from going negative
        hydraul.c:1796
        """
        pass

class Pipe(Link):
    OPEN = 1
    CLOSE = 2
    def __init__(self, pipe_data, nodes, units):

        self.length = units.length(float(pipe_data.Length))
        self.diameter = units.diam(float(pipe_data.Diameter))
        self.roughness = float(pipe_data.Roughness)

        #0.02517*Link[k].Km/SQR(Link[k].Diam)/SQR(Link[k].Diam);
        self.minor_loss = float(pipe_data.MinorLoss) # needs to be converted
        dia_sqr = self.diameter*self.diameter
        self.minor_loss = 0.02517 * self.minor_loss / dia_sqr / dia_sqr
        self.status = Pipe.OPEN
        try:
            self.status = pipe_data.Status.strip().toupper() == "CLOSED"
        except:
            pass

        Link.__init__(self, pipe_data, nodes)

        
    def initial_flow(self):
        #print "%30.16f" % (self.diameter*self.diameter*PI/4.0)
        #hydraul.c:324
        return self.diameter*self.diameter*PI/4.0

    def resistance(self):
        return 4.727*self.length*math.pow(self.roughness,-EPAnet.Hexp)*math.pow(self.diameter,-4.871) #H-W roughness coefficient

    def update(self, Q):
        q = abs(Q)
        hpipe = self.r * math.pow(q, EPAnet.Hexp)
        p = EPAnet.Hexp*hpipe
        hml = 0.0
        if self.minor_loss > 0.0:
            hml = self.minor_loss*q*q
            p += 2.0*hml

        p = Q/p
        self.id_headloss = abs(p)
        self.flow_correction_factor = p*(hpipe + hml)

    def calc_id_headloss(self, Q):
        q = abs(Q)
        p = EPAnet.Hexp * self.r * math.pow(q, EPAnet.Hexp)
        
        if self.minor_loss > 0.0:
            p += 2.0*self.minor_loss*q*q
        return abs(Q/p)

    def calc_flow_correction_factor(self, Q):
        q = abs(Q)
        y = self.r * math.pow(q, EPAnet.Hexp)
        if self.minor_loss > 0.0:
            y += 2.0*self.minor_loss*q*q
        return y*self.id_headloss(Q)

class Curve():
    def __init__(self, curve_data, id):
        self.id = id
        self.x = []
        self.y = []

        for cd in curve_data:
            if cd.ID <> self.id:
                continue
            self.x.append(float(cd.X))
            self.y.append(float(cd.Y))

class Node():
    def __init__(self, node_data):
        self.id = node_data.ID
        self.idx = None
        self.demand = 0.0
        self.initial_level = 0.0 #overwritten by tank and reservoir
        self.closed = False
        self.matrix_offset = None  #direct offset into csr data array


class Junction(Node):
    def __init__(self, junction_data, units):
        Node.__init__(self, junction_data)
        self.demand = units.from_gpm(float(junction_data.Demand))
        self.elevation = float(junction_data.Elevation) #TODO convert
        self.pattern = None
        if junction_data.Pattern is not None:
            raise NotImplemented


class Reservoir(Node):
    def __init__(self, reservoir_data, units):
        Node.__init__(self, reservoir_data)
        self.head = units.length(float(reservoir_data.Head)) #TODO convert
        self.initial_level = self.head
        self.elevation = 0.0
        self.pattern = reservoir_data.Pattern
        if self.pattern:
            raise NotImplemented



class Tank(Node):
    def __init__(self, tank_data, units):
        Node.__init__(self, tank_data)
        #Elevation   	InitLevel   	MinLevel    	MaxLevel    	Diameter    	MinVol      	VolCurve
        self.elevation = float(tank_data.Elevation)
        self.initial_level = self.elevation + float(tank_data.InitLevel)
        self.minimum_level = float(tank_data.MinLevel)
        self.maximum_level = float(tank_data.MaxLevel)


US_UNITS = [input.UNITS_CFS,
            input.UNITS_GPM,
            input.UNITS_MGD,
            input.UNITS_IMGD,
            input.UNITS_AFD]

class EPAnet():
    Hexp = 1.852 
    RQtol = 1e-7
    Hacc = 0.001 #hydraulic accuracy

    def __init__(self, input_file):
        self.__parse_input(input_file)
        self.__init_arrays()
        self.__init_matrix_offsets()

    def step(self):
        """
        perform a single calculation step
        """
        iteration = 0
        while True:
            (A, F) = self.__new_coefficients()
            self.__solve(A, F)
            relative_error = self.__new_flows()
            print "relative error: %2.14f" % relative_error
            if relative_error <= EPAnet.Hacc:
                break
                #check for maxiter
        

    def __parse_input(self, in_file):
        from input import load
        data = load(in_file)

        options = data[input.OPTIONS]

        if options.flow_unit in US_UNITS:
            self.units = UnitsUS() #TODO make proper units system
        else:
            self.units = UnitsSI() #TODO implement UnitsSI

        #init curves
        self.curves = {}
        for c in data[input.CURVES]:
            self.curves[c.ID] = None

        for id in self.curves.keys():
            self.curves[id] = Curve(data[input.CURVES], id)

        self.__prepare_nodes(data)
        self.__prepare_links(data)

    def __prepare_nodes(self, data):
        #prepare junctions
        self.junctions = collections.OrderedDict()
        for j in data[input.JUNCTIONS]:
            junction = Junction(j, self.units)
            self.junctions[junction.id] = junction

        self.reservoirs = collections.OrderedDict()
        for r in data[input.RESERVOIRS]:
            reservoir = Reservoir(r, self.units)
            self.reservoirs[reservoir.id] = reservoir

        self.tanks = collections.OrderedDict()
        for t in data[input.TANKS]:
            tank = Tank(t, self.units)
            self.tanks[tank.id] = tank

        self.nodes = collections.OrderedDict()
        self.nodes.update(self.junctions)
        self.nodes.update(self.reservoirs)
        self.nodes.update(self.tanks)

        #done with nodes here

        # set idx which is the indices into the matrix
        for i, node in enumerate(self.nodes.values()):
            node.idx = i

    def __prepare_links(self, data):
        self.pumps = collections.OrderedDict()
        for p in data[input.PUMPS]:
            pump = Pump(p, self.curves, self.nodes, self.units)
            self.pumps[pump.id] = pump


        self.pipes = collections.OrderedDict()
        for p in data[input.PIPES]:
            pipe = Pipe(p, self.nodes, self.units)
            self.pipes[pipe.id] = pipe

        self.links = self.pipes.values() + self.pumps.values() #+ valves TODO

    def __init_arrays(self):
        #init Q
        self.Q = np.array([link.initial_flow() for link in self.links])
        #print Q

        #init D
        self.D = np.array([j.demand for j in self.nodes.values()])
        #print D

        #init H
        self.H = np.array([n.initial_level for n in self.nodes.values()])
        #print H

        self.F = np.zeros(len(self.junctions))
        self.X = np.zeros(len(self.links))
        self.P = np.zeros(len(self.links))
        self.Y = np.zeros(len(self.links))

        #njuncs = len(self.junctions)
        #self.A = sparse.csc_matrix((njuncs,njuncs))

    def __init_matrix_offsets(self):
        njuncs = len(self.junctions)
        matrix = sparse.lil_matrix((njuncs,njuncs))
        for link in self.links:
            n1 = link.node1
            n2 = link.node2

            n1junc = isinstance(n1, Junction)
            n2junc = isinstance(n2, Junction)

            if n1junc and n2junc:
                matrix[n1.idx, n2.idx] = 1.0
                matrix[n2.idx, n1.idx] = 1.0
            if n1junc:
                matrix[n1.idx, n1.idx] = 1.0
            if n2junc:
                matrix[n2.idx, n2.idx] = 1.0

        self.A = matrix.tocsr()

        for link in self.links:
            n1 = link.node1
            n2 = link.node2

            n1junc = isinstance(n1, Junction)
            n2junc = isinstance(n2, Junction)

            if n1junc and n2junc:
                upper_found, lower_found = (False, False)
                for c in range(self.A.indptr[n1.idx], self.A.indptr[n1.idx+1]):
                    if self.A.indices[c] == n2.idx:
                        upper_found = True
                        link.upper_matrix_offset = c

                for c in range(self.A.indptr[n2.idx], self.A.indptr[n2.idx+1]):
                    if self.A.indices[c] == n1.idx:
                        lower_found = True
                        link.lower_matrix_offset = c

                assert(upper_found and lower_found)

            if n1junc:
                found = False
                for c in range(self.A.indptr[n1.idx], self.A.indptr[n1.idx+1]):
                    if self.A.indices[c] == n1.idx:
                        n1.matrix_offset = c
                        found = True
                assert found

            if n2junc:
                found = False
                for c in range(self.A.indptr[n2.idx], self.A.indptr[n2.idx+1]):
                    if self.A.indices[c] == n2.idx:
                        n2.matrix_offset = c
                        found = True
                assert found
            self.A.data.fill(0.0)

    def __new_coefficients(self):
        n = len(self.nodes.values())
        njuncs = len(self.junctions)

        self.F.fill(0.0)
        self.X.fill(0.0)
        self.P.fill(0.0)
        self.Y.fill(0.0)
        self.A.data.fill(0.0)

        self.__new_link_coefficients(self.A, self.X, self.F)
        self.__new_emitter_coefficients()
        self.__new_node_coefficients(self.X, self.F)
        self.__new_valve_coefficients()

        return (self.A, self.F)


    def __new_link_coefficients(self, A, X, F):
        for i, link in enumerate(self.links):
            n1 = link.node1
            n2 = link.node2

            #if not isinstance(n1, Junction) or not isinstance(n2, Junction):
            #    print "nojunccon: ", (n1.idx, n2.idx), n1, n2
            #if isinstance(n1, Reservoir) or isinstance(n2, Reservoir):
            #    print "Tankcon:", (n1.idx, n2.idx)

            X[n1.idx] -= self.Q[i]
            X[n2.idx] += self.Q[i]

            link.update(self.Q[i]) #calc id_headloss and flow_correction_factor
            self.P[i] = link.id_headloss
            self.Y[i] = link.flow_correction_factor

            n1junc = isinstance(n1, Junction)
            n2junc = isinstance(n2, Junction)

            if n1junc and n2junc:
                A.data[link.upper_matrix_offset] -= self.P[i]
                A.data[link.lower_matrix_offset] -= self.P[i]

            if n1junc:
                A.data[n1.matrix_offset] += self.P[i]
                F[n1.idx] += self.Y[i]
            else:
                F[n2.idx] += self.P[i]*self.H[n1.idx]


            if n2junc:
                A.data[n2.matrix_offset] += self.P[i]
                F[n2.idx] -= self.Y[i]
            else:
                F[n1.idx] += self.P[i]*self.H[n2.idx]


    def __new_emitter_coefficients(self):
        #TODO implement
        pass

    def __new_valve_coefficients(self):
        #TODO implement
        pass

    def __new_node_coefficients(self, X, F):
        for i in xrange(len(self.junctions)):
            X[i] -= self.D[i]
            F[i] += X[i]

    def __solve(self, A, F):
        njuncs = len(self.junctions)
        self.H[0:njuncs] = splinalg.spsolve(A, F[:njuncs])

    def __new_flows(self):

        for i, node in enumerate(self.nodes.values()):
            #print node
            if isinstance(node, Tank) or isinstance(node, Reservoir):
                self.D[i] = 0.0

        qsum = 0.0
        dqsum = 0.0

        for i, link in enumerate(self.links):
            n1 = link.node1
            n2 = link.node2

            dh = self.H[n1.idx] - self.H[n2.idx]
            dq = self.Y[i] - self.P[i]*dh

            if isinstance(link, Pump):
                link.correct_flow()

            self.Q[i] -= dq

            qsum += abs(self.Q[i])
            dqsum += abs(dq)

            if not n1.closed and (isinstance(n1, Tank) or isinstance(n1, Reservoir)):
                #TODO check not closed
                self.D[n1.idx] -= self.Q[i]

            if not n2.closed and (isinstance(n2, Tank) or isinstance(n2, Reservoir)):
                #TODO check not closed
                self.D[n2.idx] += self.Q[i]

        #TODO update emitter flows

        if qsum > EPAnet.Hacc:
            return dqsum/qsum

        return dqsum




def main(argv):
    epanet = EPAnet(argv[1])
    epanet.step()

if __name__ == "__main__":
    assert(len(sys.argv) > 1)
    profile_file = "/tmp/epanetprof.txt"
    import cProfile
    cProfile.run("main(sys.argv)", profile_file)
    import pstats
    #stats = pstats.Stats(profile_file)
    #stats.sort_stats('time', ).print_stats()
