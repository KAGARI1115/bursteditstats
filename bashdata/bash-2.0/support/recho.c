#include <stdio.h>

main(argc, argv)
int	argc;
char	**argv;
{
	register int	i;

	for (i = 1; i < argc; i++) {
		printf("argv[%d] = <", i);
		strprint(argv[i]);
		printf(">\n");
	}
}

strprint(str)
char	*str;
{
	register char *s;
	int	c;

	for (s = str; s && *s; s++) {
		if (*s < ' ') {
			putchar('^');
			putchar(*s+64);
		} else if (*s == 127) {
			putchar('^');
			putchar('?');
		} else
			putchar(*s);
	}
	return(0);
}
