// note: originally written by Robert Koenighofer.

#ifndef Logger_H__
#define Logger_H__


#include <iostream>
#include <ostream>
#include <vector>
#include <cstdlib>


using namespace std;


class Logger {
public:
    static Logger &instance();

    enum Logtype {
        ERR = 0,
        WRN = 1,
        RES = 2,
        INF = 3,
        DBG = 4,
        LOG = 5
    };

    void enable(Logtype lt);
    void disable(Logtype lt);
    bool isEnabled(Logtype lt) const;

    ostream &getOstream(Logtype lt);

protected:
    static Logger *instance_;
    vector<ostream *> streams_;
    vector<bool> enabled_;

private:
    Logger();
    ~Logger();

    Logger(const Logger &other);
    Logger &operator=(const Logger &other);
};

#define L_ERR(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::ERR))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::ERR)) <<                       \
    "[ERR] " << message << std::endl;                                     \
  }                                                                       \
}

#define L_WRN(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::WRN))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::WRN)) <<                       \
    "[WRN] " << message << std::endl;                                     \
  }                                                                       \
}

#define L_RES(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::RES))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::RES)) <<                       \
    "[RES] " << message << std::endl;                                     \
  }                                                                       \
}

#define L_INF(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::INF))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::INF)) <<                       \
    "[INF] " << message << std::endl;                                     \
  }                                                                       \
}

#define L_DBG(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::DBG))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::DBG)) <<                       \
    "[DBG] " << message << std::endl;                                     \
  }                                                                       \
}

#define L_LOG(message)                                                    \
{                                                                         \
  if (Logger::instance().isEnabled(Logger::LOG))                          \
  {                                                                       \
    (Logger::instance().getOstream(Logger::LOG)) <<                       \
    "[LOG] " << message << std::endl;                                     \
  }                                                                       \
}


#define MASSERT(condition, message)                                          \
{                                                                            \
  if(!(condition))                                                           \
  {                                                                          \
    std::cerr << __FILE__ << " (" << __LINE__ << ") : " << message << endl;  \
    abort();                                                                 \
  }                                                                          \
}


#endif // Logger_H__
