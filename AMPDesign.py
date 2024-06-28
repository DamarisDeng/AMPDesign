import getopt
import sys
from colorama import Fore
from models.seqgan.Seqgan import Seqgan

def set_gan(output_path):
    gan = Seqgan(output_path)
    gan.vocab_size = 5000
    gan.generate_num = 10000
    return gan

def parse_cmd(argv):
    usr_msg = 'usage: python AMPDesign.py -s <species> -a <model-path> -d train.fasta -o output.fasta'
    try:
        opts, args = getopt.getopt(argv, "hs:a:d:o:")
        opt_arg = dict(opts)
        if '-h' in opt_arg:
            print(usr_msg)
            sys.exit(0)
        required_opts = ['-a', '-d', '-o']
        for opt in required_opts:
            if opt not in opt_arg:
                print(f"Missing required argument: {opt}")
                print(usr_msg)
                sys.exit(1)
        gan = set_gan(opt_arg['-o'])  # specified output file
        gan.train_real(data_loc=opt_arg['-d'], model_loc=opt_arg['-a'], output_path=opt_arg['-o'])


    except getopt.GetoptError:
        print('Invalid arguments!')
        print('`python AMPDesign.py -h` for help')
        sys.exit(-1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(-1)

if __name__ == '__main__':
    parse_cmd(sys.argv[1:])
