#include "Logger.hpp"


Logger *Logger::instance_ = NULL;

Logger &Logger::instance()
{
  if(instance_ == NULL)
    instance_ = new Logger;
  MASSERT(instance_ != NULL, "could not create Logger instance");
  return *instance_;
}

Logger::Logger ()
{
  streams_.push_back(&cerr); // ERR
  streams_.push_back(&cerr); // WRN
  streams_.push_back(&cout); // RES
  streams_.push_back(&cout); // INF
  streams_.push_back(&cout); // DBG
  streams_.push_back(&cout); // LOG
  enabled_.push_back(true);  // ERR
  enabled_.push_back(true);  // WRN
  enabled_.push_back(true);  // RES
  enabled_.push_back(true);  // INF
  enabled_.push_back(false); // DBG
  enabled_.push_back(true);  // LOG
}

Logger::~Logger ()
{
  streams_[0]->flush();
  streams_[1]->flush();
  streams_[2]->flush();
  streams_[3]->flush();
  streams_[4]->flush();
  streams_[5]->flush();
}

void Logger::disable(Logtype lt)
{
  enabled_[lt] = false;
}

void Logger::enable(Logtype lt)
{
  enabled_[lt] = true;
}

bool Logger::isEnabled(Logtype lt) const
{
  return enabled_[lt];
}

ostream &Logger::getOstream(Logtype lt)
{
  MASSERT(streams_[lt]->good(), "Cannot write to ostream");
  return *streams_[lt];
}

