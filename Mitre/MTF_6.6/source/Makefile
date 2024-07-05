#                                 NOTICE
#                                 ------
# 
#                This "sinemtf" software was developed by 
#                         The MITRE Corporation
# and was produced for the U.S. Government under Contract numbers 
# J-FBI-12-128, J-FBI-86-028, J-FBI-93-039, DAAB07-99-C-C201, DAAB07-00-C-C201, 
# DAAB07-01-C-C201, DAAB07-01-C-N200, DAAB07-02-C-N200, W15P7T-04-C-D001,
# and W15P7T-04-C-D199 and is subject to the Rights in Data-General 
# clause FAR 52.227-l4, Alt. IV (DEC 2007). 
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# -Redistribution of source code must retain the above copyright notice, this
# list of conditions, and the following disclaimer.
# -Redistribution in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

# To compile a version of the mtf program with the TIFF software:
# - make sure that -DUSE_TIFF is in the CPPFLAGS and
# - make sure that -ltiff in LDLIBS
# - as well as any other directives to find TIFF includes and libraries.

# To compile a version of the mtf program without the TIFF software:
# - remove -DUSE_TIFF from the CPPFLAGS
# - remove -ltiff from LDLIBS

# On Windows machines, make sure to define MSDOS (CPPFLAGS += -DMSDOS)

# Set the appropriate compiler and basic switches
CC = gcc
CFLAGS = -O3 -pedantic -ansi -Wall
#CFLAGS = -g -pedantic -ansi -Wall

# Enable TIFF usage 
LDLIBS = -ltiff 
#LDLIBS =
CPPFLAGS = -DUSE_TIFF 

# For Windows Systems: (Others comment out following line)
CPPFLAGS += -DMSDOS

OBJS = \
	density.o \
	mtf.o \
	scurve.o \
	sine.o \
	spline.o \
	skew.o \
	corner.o


${OBJS}: mtf.h xterndef.h xterndcl.h 

default: mtf

mtf:	$(OBJS)
	$(LINK.c) -o $@ $(OBJS) $(LDLIBS) -lm

clean: 
	rm *.o core