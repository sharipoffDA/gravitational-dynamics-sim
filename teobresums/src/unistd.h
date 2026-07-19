#ifndef UNISTD_H
#define UNISTD_H

/* Этот файл — заглушка POSIX для MSVC */
#ifdef _WIN32
    #include <io.h>
    #include <process.h>
    #include <direct.h>

    /* Некоторые стандартные POSIX константы, если они понадобятся библиотеке */
    #define F_OK 0
    #define W_OK 2
    #define R_OK 4
#endif

#endif /* UNISTD_H */