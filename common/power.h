/*
 power.h

 This is part of the tree library

 Copyright 2015 Ibrahim Umar (UiT the Arctic University of Norway)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */


#ifndef power_h
#define power_h


#if (defined (__x86_64__) && !defined(__KNC__))

#include "pcmpower.h"

#define ENERGY_START() pcm_bench_start();
#define ENERGY_END()  {pcm_bench_end(); pcm_bench_print();}

#endif


#ifdef __arm__

#include "armpower.h"

#define ENERGY_START()  armpower_start();
#define ENERGY_END()  armpower_finalize();

#endif

#ifdef __KNC__

#include "micpower.h"

#define ENERGY_START()  micpower_start();
#define ENERGY_END()  micpower_finalize();

#endif

#endif
