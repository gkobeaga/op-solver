int
prime_check(unsigned int p)
{
    unsigned int i;

    if ((p & 1) == 0)
        return 0;
    for (i = 3; i * i <= p; i += 2)
    {
        if (p % i == 0)
            return 0;
    }
    return 1;
}

unsigned int
prime_next(unsigned int x)
{
    if (x < 3)
        return 3;
    x |= 1;
    while (!prime_check(x)) x += 2;
    return x;
}
