target_sources(PhysicaLogger PUBLIC FILE_SET CXX_MODULES FILES
    LoggerImpl/AbstractLogger.cppm
    LoggerImpl/FileLogger.cppm
    LoggerImpl/LogBuffer.cppm
    LoggerImpl/LoggerTimer.cppm
    LoggerImpl/StdLogger.cppm
    FormatAnalyzer.cppm
    Logger.cppm
    LoggerType.cppm
)
