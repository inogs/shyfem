#!/bin/sh

#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

gfortran nudge.f
[ $? -ne 0 ] && exit 1

./a.out

