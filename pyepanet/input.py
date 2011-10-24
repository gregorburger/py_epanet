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

import sys
import os

HLOSS_CM = "C-M"

HLOSS_DW = "D-W"

HLOSS_HW = "H-W"

OPTIONS = 'OPTIONS'

UNITS_CMD = "CMD"
UNITS_CMH = "CMH"
UNITS_MLD = "MLD"
UNITS_LPM = "LPM"
UNITS_LPS = "LPS"
UNITS_AFD = "AFD"
UNITS_IMGD = "IMGD"
UNITS_MGD = "MGD"
UNITS_GPM = "GPM"
UNITS_CFS = "CFS"

RESERVOIRS = 'RESERVOIRS'
PUMPS = 'PUMPS'
PIPES = 'PIPES'
TANKS = 'TANKS'
JUNCTIONS = 'JUNCTIONS'
CURVES = 'CURVES'


def is_empty(line):
    return len(line.strip()) == 0

def is_section(line):
    return line.strip()[0] == '[' and line.strip()[-1] == ']'

def process_title(section, line, data):
    if section not in data:
        data[section] = line
    else:
        data[section] += line

class SectionEntry:
    def __init__(self, entries, line):
        values = line.strip().split()
        #print "values %s" % len(values)
        #print "entries %s" % len(entries)
        #print ""
        if len(values) < len(entries):
            values += [None] * (len(entries) - len(values)) # append optionals as None values
        assert len(values) == len(entries)
        self.__dict__.update(dict(zip(entries, values)))

class KVSectionEntry:
    def __init__(self, key, value):
        self.__dict__[key] = value


def process_simple_section(section, line, data, entries):
    se = SectionEntry(entries, line)
    if section not in data:
        data[section] = [se]
    else:
        data[section].append(se)

def process_kv_section(section, line, data):
    kv = line.strip().split(None, 1)
    print kv
    assert (len(kv) == 2)
    se = KVSectionEntry(kv[0], kv[1])

    if section not in data:
        data[section] = [se]
    else:
        data[section].append(se)

    print "kv section %s" % section

def process_demands(section, line, data):
    pass


def process_times(section, line, data):
    pass

class PumpSectionEntry:
    def __init__(self, line):
        items = line.strip().split()
        kw = ['POWER', 'HEAD', 'PATTERN', 'SPEED']
        self.ID = items[0]
        self.Node1 = items[1]
        self.Node2 = items[2]
        try:
            value = float(items[3])
            print value
            #print "version 1 pump section"
        except ValueError:
            #print "version 2 pump section"
            parameters = {}
            if items[3] not in kw:
                print "unknown pumps section format"
            for i in xrange(3, 3 + (len(items)-3)/2):
                key, value = items[i], items[i+1]
                #print key, value
                parameters[key] = value

            self.parameters = parameters

def process_pumps(section, line, data):
    se = PumpSectionEntry(line)
    if section not in data:
        data[section] = [se]
    else:
        data[section].append(se)




class Options:
    def __init__(self):
        self.flow_unit = UNITS_GPM
        self.headloss_method = HLOSS_HW
        self.hydraulics_save_file = ""
        self.hydraulics_use_file = ""

UNITS = [UNITS_CFS,
         UNITS_GPM,
         UNITS_MGD,
         UNITS_IMGD,
         UNITS_AFD,
         UNITS_LPS,
         UNITS_LPM,
         UNITS_MLD,
         UNITS_CMH,
         UNITS_CMD]

HEADLOSS_METHODS = [HLOSS_HW,
                    HLOSS_DW,
                    HLOSS_CM]

def capitalize(str):
    str.capitalize()

def process_options(section, line, data):
    if section not in data:
        data[section] = Options()

    items = line.strip().split()
    
    cap_items = line.strip().split()
    map(capitalize, cap_items)

    if cap_items[0] == "UNITS" and cap_items[1] in UNITS:
        options.flow_unit = cap_items[1]

    if cap_items[0] == "HEADLOSS" and cap_items[1] in HEADLOSS_METHODS:
        options.headloss_method = cap_items[1]

    if cap_items[0] == "HYDRAULICS":
        assert(cap_items[1] in ["USE", "SAVE"])
        if cap_items[1] == "USE":
            options.hydraulics_use_file = items[2]
        if cap_items[1] == "SAVE":
            options.hydraulics_save_file = items[2]


SECTIONS = {
    'TITLE' : process_title,
    #'VALVES' : process_valves,
    'DEMANDS' : process_demands,
    'TIMES' : process_times,
    OPTIONS: process_options,
    'PUMPS' : process_pumps,
    }

SIMPLE_SECTIONS = {
    JUNCTIONS: ["ID", "Elevation", "Demand", "Pattern"],
    RESERVOIRS: ["ID", "Head", "Pattern"],
    TANKS: ["ID", "Elevation", "InitLevel", "MinLevel", "MaxLevel", "Diameter", "MinVol", "VolCurve"],
    PIPES: ["ID", "Node1", "Node2", "Length", "Diameter", "Roughness", "MinorLoss", "Status"],
    'VALVES' : ["ID", "Node1", "Node2", "Diameter", "Type", "Setting", "MinorLoss"],
    CURVES: ['ID', 'X', 'Y']
    #PUMPS: ["ID", "Node1", "Node2", "Pcurve", "Status"],
}

IGNORED_SECTIONS = {
    'PATTERNS',
    'CURVES',
    'CONTROLS',
    'ENERGY',
    'QUALITY',
    'REACTIONS',
    'COORDINATES',
    'LABELS',
    'BACKDROP',
    'REPORT'
}

KV_SECTIONS = {
    #'OPTIONS'
}

def process_for_section(section, line, data):
    if section in SIMPLE_SECTIONS:
        process_simple_section(section, line, data, SIMPLE_SECTIONS[section])
        return
    if section in KV_SECTIONS:
        process_kv_section(section, line, data)
        return
    if section in SECTIONS:
        SECTIONS[section](section, line, data)
        return
    if section in IGNORED_SECTIONS:
        return

    print 'section %s not implemented' % section
    exit()

def load(path):
    data = {}
    with open(path) as f:
        current_section = ''


        for line in f:
            #strip comments
            line = line.partition(';')[0]

            if is_empty(line):
                continue

            if is_section(line):
                current_section = line.strip()[1:-1]
                #print "found section: %s" % current_section
                continue

            if current_section:
                #print 'in_section:'
                process_for_section(current_section, line, data)
                continue

            print 'unknown line %s' % line.strip()

    return data



if __name__ == "__main__":
    data = load(sys.argv[1])
    print 'title of the project is: %s' % data['TITLE']

    for j in data[JUNCTIONS]:
        print "%s %s %s" % (j.ID, j.Elevation, j.Demand)

