/* From Loic Grenie */
/* to build
gcc -o cook -O2 cook.c
cp cook uncook
to use (see mpigp)
./uncook mpirun -np 1 ./cook ./gp : -np 4 ./gp
*/

#define _XOPEN_SOURCE 600
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/select.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <termios.h>
#include <sys/ioctl.h>

int dead = 0, status;

void chld_handler(int sig)
{
    waitpid(-1, &status, WUNTRACED);
    if (WIFEXITED(status))
        dead = 1;
    else if (WIFSIGNALED(status))
        dead = 1;
    else if (WIFSTOPPED(status))
        kill(getpid(), WSTOPSIG(status));
}

void usage(char *argv0, int exitok)
{
    fprintf(exitok ? stdout : stderr,
            "Usage: %s [-u] [-c] [--] program [args]\n", argv0);
    exit(!exitok);
}

int main(int argc, char **argv)
{
    int cook, tiook = 0;
    struct termios tio;
    int fd;
    char *name, *argv0 = argv[0];
    char bufin[1024], bufout[1024];
    int lgin = 0, lgout = 0;

    if (argv[0]) {
        char *p = strrchr(argv[0], '/');

        if (!p)
            p = argv[0];
        else
            p++;
        cook = strcmp(p, "uncook");
    }
    argc--;
    argv++;
    while (argc > 0 && *argv && **argv == '-') {
        char *p = *argv;
        if (!strcmp(*argv, "--")) {
            argc--;
            argv++;
            break;
        }
        else if (!strcmp(*argv, "--help"))
            usage(argv0, 1);
        while (*++p) {
            switch (*p) {
            case 'u':
                cook = 0;
                break;
            case 'c':
                cook = 1;
                break;
            default:
                usage(argv0, *p == 'h');
                break;
            }
        }
    }
    if (argc <= 0 || !argv)
        usage(argv0, 0);
    if (cook) {
        fd = posix_openpt(O_RDWR);
        if (fd == 0) {
            fd = dup(fd);
            if (fd == 1) {
                fd = dup(fd);
                close(1);
            }
            close(0);
        }
        if (grantpt(fd) < 0) {
            perror("grantpt");
            exit(3);
        }
        if (unlockpt(fd)) {
            perror("unlockpt");
            exit(4);
        }
        name = ptsname(fd);
    }
    else {
        struct winsize wsz;

        if (!ioctl(0, TIOCGWINSZ, &wsz)) {
            char buf[1024];
            sprintf(buf, "%d", wsz.ws_row);
            setenv("ROWS", buf, 0);
            sprintf(buf, "%d", wsz.ws_col);
            setenv("COLUMNS", buf, 0);
        }
        if (!tcgetattr(0, &tio))
        {
            struct termios tio2;

            tiook = 1;
            bcopy(&tio, &tio2, sizeof tio);
            tio2.c_iflag &= ~(IGNBRK | BRKINT | PARMRK | ISTRIP
                           | INLCR | IGNCR | ICRNL | IXON);
            tio2.c_oflag &= ~OPOST;
            tio2.c_lflag &= ~(ECHO | ECHONL | ICANON | ISIG | IEXTEN);
            tio2.c_cflag &= ~(CSIZE | PARENB);
            tio2.c_cflag |= CS8;
            tcsetattr(0, TCSADRAIN, &tio2);
        }
        else
            (void)system("stty -isig -icanon -echo");
    }
    switch(fork()) {
    case -1:
        perror("fork");
        exit(1);
    case 0:
        /* Child */
        if (cook) {
            char *prows, *pcols;
            struct winsize wsz;

            close(0); close(1); close(2);
            setsid();
            close(fd);
            open(name, O_RDWR); /* stdin */
            open(name, O_RDWR); /* stdout */
            (void)system("stty sane pass8"); /* Before opening stderr */
            open(name, O_RDWR); /* stderr */
            if ((prows = getenv("ROWS")) && (pcols = getenv("COLUMNS"))) {
                wsz.ws_row = atoi(prows);
                wsz.ws_col = atoi(pcols);
                ioctl(0, TIOCSWINSZ, &wsz);
            }
        }
        execvp(argv[0], argv);
        fprintf(stderr, "Command not found\n");
        exit(2);
    }
    signal(SIGCHLD, chld_handler);

    if (cook) {
        fcntl(0, F_SETFL, fcntl(0, F_GETFL, 0) & ~O_NONBLOCK);
        fcntl(fd, F_SETFL, fcntl(fd, F_GETFL, 0) & ~O_NONBLOCK);
        while (!dead || lgout) {
            int maxfd = 2, r;
            fd_set inset, outset;

            FD_ZERO(&inset);
            FD_ZERO(&outset);
            if (lgin < sizeof(bufin)) {
                FD_SET(0, &inset);
                if (maxfd < 1) maxfd = 1;
            }
            if (lgout) {
                FD_SET(1, &outset);
                if (maxfd < 2) maxfd = 2;
            }
            if (!dead && lgin) {
                FD_SET(fd, &outset);
                if (maxfd <= fd) maxfd = fd + 1;
            }
            if (!dead && lgout < sizeof(bufout)) {
                FD_SET(fd, &inset);
                if (maxfd <= fd) maxfd = fd + 1;
            }
            select(maxfd, &inset, &outset, NULL, NULL);
            if (FD_ISSET(0, &inset)) {
                r = read(0, bufin+lgin, sizeof(bufin) - lgin);
                if (r > 0)
                    lgin += r;
                else
                    close(0);
            }
            if (FD_ISSET(1, &outset)) {
                r = write(1, bufout, lgout);
                if (r <= 0)
                    lgout = 0;
                else
                    lgout -= r;
            }
            if (FD_ISSET(fd, &inset)) {
                r = read(fd, bufout+lgout, sizeof(bufout) - lgout);
                if (r <= 0) exit(0);
                lgout += r;
            }
            if (FD_ISSET(fd, &outset)) {
                r = write(fd, bufin, lgin);
                if (r <= 0) exit(0);
                lgin -= r;
            }
        }
    }
    else {
        while (!dead)
            chld_handler(0);
    }
    if (tiook)
        tcsetattr(0, TCSADRAIN, &tio);
    if (dead) {
        if (WIFEXITED(status))
            exit(WEXITSTATUS(status));
        else if (WIFSIGNALED(status))
            kill(getpid(), WTERMSIG(status));
        fprintf(stderr, "Unknown exit status %08x\n", status);
        exit(1);
    }
    /* NOT REACHED */
    exit(0);
}
