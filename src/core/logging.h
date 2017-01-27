#pragma once

// Homegrown logging
//

#include <iostream>

#define LOGGING_CERR if(true) {} else std::cerr

#define CHECK(a) if(!(a)) LOGGING_CERR
#define CHECK_EQ(a, b)
#define CHECK_GE(a, b)
#define CHECK_GT(a, b)
#define CHECK_LE(a, b) if(!((a) <= (b))) LOGGING_CERR
#define CHECK_LT(a, b) if(!((a) < (b))) LOGGING_CERR
#define CHECK_NE(a, b)
#define CHECK_NOTNULL(a)
#define DCHECK(a)
#define DCHECK_LT(a,b)
#define DCHECK_NE(a,b)
#define FATAL
#define LOG(a) LOGGING_CERR
#define VLOG(a) LOGGING_CERR

extern bool FLAGS_logtostderr;
extern int FLAGS_minloglevel;
extern int FLAGS_stderrthreshold;
extern int FLAGS_v;
extern char* FLAGS_log_dir;
