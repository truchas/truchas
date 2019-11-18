/*******************************************************************************
* Copyright Notice
*  + 2010-2012 North Carolina State University
*  + 2010-2012 Pacific Northwest National Laboratory
* 
* This file is part of SCORPIO.
* 
* SCORPIO is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
* 
* SCORPIO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with SCORPIO.  If not, see <http://www.gnu.org/licenses/>.
* 
*******************************************************************************/
/** 
 * @file scorpiof.h 
 * @brief Fortran Header file for the parallel I/O framework interface.
 * @author Sarat Sreepathi (sarat@computer.org)
 * @version 0.01
 * @date 2010-11-01
 */

#ifndef _PARALLELIOF_H
#define _PARALLELIOF_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "scorpio.h"

// For now, just define this flag here
// #define LINUX_FTN
// As we get more platforms, we can expand this section 
// #define IBM_FTN
#if defined(IBM_FTN)
#define F2C(name) f ## name 
#elif defined(LINUX_FTN)
// Concatenate _ to end of name
#define F2C(name) f ## name ## _
#endif


	// typedef iogroup_t F2C(iogroup_t); 
	// typedef file_mode_t F2C(file_mode_t);
	// typedef iopattern_t F2C(iopattern_t);
	// typedef datatype_t F2C(datatype_t);

	extern iogroup_t **iogroups;



#ifdef __cplusplus
}
#endif

#endif

