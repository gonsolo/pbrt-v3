#pragma once

#include <iostream>

#define LOGGING_CERR if(true) {} else std::cerr

class DoNothing {};
static DoNothing doNothing;

template<typename T>
DoNothing CHECK(const T& t) { return doNothing; }

template<typename T>
DoNothing CHECK_NOTNULL(const T& t) { return doNothing; }

template<typename T>
DoNothing DCHECK(const T& t) { return doNothing; }

template<typename T>
DoNothing LOG(const T& t) { return doNothing; }

template<typename T>
DoNothing VLOG(const T& t) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_EQ(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_GE(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_GT(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_LE(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_LT(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing CHECK_NE(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing DCHECK_LT(const T& a, const U& b) { return doNothing; }

template<typename T, typename U>
DoNothing DCHECK_NE(const T& a, const U& b) { return doNothing; }


template<typename T>
inline DoNothing operator<<(DoNothing dn, T value) { return dn; }


extern bool  FLAGS_logtostderr;
extern int   FLAGS_minloglevel;
extern int   FLAGS_stderrthreshold;
extern int   FLAGS_v;
extern char* FLAGS_log_dir;

const int INFO    = 0;
const int WARNING = 1;
const int ERROR   = 2;
const int FATAL   = 3;

