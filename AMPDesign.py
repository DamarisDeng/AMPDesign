import getopt
import sys
from colorama import Fore
from models.seqgan.Seqgan import Seqgan


def set_gan(gan_name):
    gans = dict()
    gans['seqgan'] = Seqgan
    Gan = gans[gan_name.lower()]
    gan = Gan()
    gan.vocab_size = 5000
    gan.generate_num = 10000
    return gan


def set_training(gan, training_method):
    try:
        gan_func = gan.train_real
    except AttributeError:
        print(Fore.RED + 'Unsupported training setting: ' + training_method + Fore.RESET)
        sys.exit(-3)
    return gan_func


def parse_cmd(argv):
    try:
        opts, args = getopt.getopt(argv, "hs:a:d:o:")
        opt_arg = dict(opts)
        if '-h' in opt_arg.keys():
            print('usage: python3 AMPDesign.py -s <microbial-type> -a <model-location>  -d train.fasta -o output.fasta')
            sys.exit(0)
        set_gan('seqgan')
        gan_func = set_training(gan, 'real') 
        if '-d' in opt_arg.keys():
            gan_func(opt_arg['-d']) # should add -s, -a, -o
            
    except getopt.GetoptError:
        print('invalid arguments!')
        print('`python AMPDesign.py -h`  for help')
        sys.exit(-1)
    pass


if __name__ == '__main__':
    gan = 'seqgan'
    parse_cmd(sys.argv[1:])
