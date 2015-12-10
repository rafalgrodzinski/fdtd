#include "utils.h"

#include <stdlib.h>
#include <stdio.h>


void check(int condition, char *message)
{
    if (condition) {
        printf("Error: %s", message);
        exit(EXIT_FAILURE);
    }
}
