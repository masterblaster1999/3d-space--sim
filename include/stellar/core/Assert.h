#pragma once

#include <cstdlib>
#include <iostream>

#ifndef STELLAR_ASSERT
  #define STELLAR_ASSERT(expr)     do {       if (!(expr)) {         std::cerr << "STELLAR_ASSERT failed: " #expr "\n"                   << "  at " << __FILE__ << ":" << __LINE__ << "\n";         std::abort();       }     } while (0)
#endif

#ifndef STELLAR_ASSERT_MSG
  #define STELLAR_ASSERT_MSG(expr, msg)     do {       if (!(expr)) {         std::cerr << "STELLAR_ASSERT failed: " #expr "\n"                   << "  message: " << (msg) << "\n"                   << "  at " << __FILE__ << ":" << __LINE__ << "\n";         std::abort();       }     } while (0)
#endif
