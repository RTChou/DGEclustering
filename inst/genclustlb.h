/* GenClust
 * Copyright (C) 2002,2004 
 *  Vito Di Gesù, Raffaele Giancarlo,
 *  Giosuè Lo Bosco, Alessandra Raimondi
 *  and Davide Scaturro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "time.h"

#define INITFILE "out.tmp"
#define TRUE 1
#define FALSE 0


typedef unsigned char ubyte;
typedef unsigned int uint;

//type used for managing the centroids
typedef double *bartype; 

//ptype is the structure which stores the value of the centroid,
//the internal variance, the number of items
typedef struct
{
 bartype bar; 
 bartype var;
 unsigned int np;
}ptype;

typedef ptype *partitiontype;

ubyte *genclust(uint nr, uint nfeat, uint nstepmax, uint nparti, double *v, double pcc, double pmm, ubyte *pinit, ubyte bestout);
