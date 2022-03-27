#include <stdio.h> /* Parâmetros */
#define f(x) ((x) * (x))
#define numberRects 50
#define lowerLimit 2.0
#define upperLimit 5.0

int main(int argc, char *argv[])
{
    int i;
    double area, at, height, width;
    area = 0.0;
    width = (upperLimit - lowerLimit) / numberRects;
    for (i = 0; i < numberRects; i++)
    {
        at = lowerLimit + i * width + width / 2.0;
        height = f(at);
        area = area + width * height;
    }
    printf("A area entre %f e %f e: %f\n", lowerLimit, upperLimit, area);
    return 0;
}