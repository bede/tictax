import argh

from tictax import tictax


def main():
	tictax.parser.dispatch()

if __name__ == '__main__':
	main()