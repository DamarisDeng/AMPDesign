import json
import os
import sys
from time import time

from models.Gan import Gan
from models.seqgan.SeqganDataLoader import DataLoader, DisDataloader
from models.seqgan.SeqganDiscriminator import Discriminator
from models.seqgan.SeqganGenerator import Generator
from models.seqgan.SeqganReward import Reward
from utils.metrics.Nll import Nll
from utils.oracle.OracleLstm import OracleLstm
from utils.utils import *


def transform_file(input_file):
    input_dir = os.path.dirname(input_file)
    output_file = os.path.join(input_dir, 'tmp.txt')
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            transformed_line = line.strip().lower()
            transformed_line = ' '.join(transformed_line)
            outfile.write(transformed_line + '\n')


class Seqgan(Gan):
    def __init__(self, output_path: str):
        super().__init__()
        # you can change parameters, generator here
        self.vocab_size = 20
        self.emb_dim = 32
        self.hidden_dim = 32
        self.sequence_length = 20
        self.filter_size = [2, 3]
        self.num_filters = [100, 200]
        self.l2_reg_lambda = 0.2
        self.dropout_keep_prob = 0.75
        self.batch_size = 64
        self.generate_num = 128
        self.start_token = 0

        self.oracle_file = 'save/oracle.txt'
        self.generator_file = 'save/generator.txt'
        self.output_path = output_path

    def init_metric(self):
        nll = Nll(data_loader=self.oracle_data_loader, rnn=self.oracle, sess=self.sess)
        self.add_metric(nll)

        inll = Nll(data_loader=self.gen_data_loader, rnn=self.generator, sess=self.sess)
        inll.set_name('nll-test')
        self.add_metric(inll)

        from utils.metrics.DocEmbSim import DocEmbSim
        docsim = DocEmbSim(oracle_file=self.oracle_file, generator_file=self.generator_file,
                           num_vocabulary=self.vocab_size)
        self.add_metric(docsim)

    def train_discriminator(self):
        generate_samples(self.sess, self.generator, self.batch_size, self.generate_num, self.generator_file)
        self.dis_data_loader.load_train_data(self.oracle_file, self.generator_file)
        for _ in range(3):
            self.dis_data_loader.next_batch()
            x_batch, y_batch = self.dis_data_loader.next_batch()
            feed = {
                self.discriminator.input_x: x_batch,
                self.discriminator.input_y: y_batch,
            }
            loss, _ = self.sess.run([self.discriminator.d_loss, self.discriminator.train_op], feed)
            print(loss)

    def evaluate(self):
        generate_samples(self.sess, self.generator, self.batch_size, self.generate_num, self.generator_file)
        if self.oracle_data_loader is not None:
            self.oracle_data_loader.create_batches(self.generator_file)
        if self.log is not None:
            if self.epoch == 0 or self.epoch == 1:
                self.log.write('epochs, ')
                for metric in self.metrics:
                    self.log.write(metric.get_name() + ',')
                self.log.write('\n')
            scores = super().evaluate()
            for score in scores:
                self.log.write(str(score) + ',')
            self.log.write('\n')
            return scores
        return super().evaluate()

    def init_oracle_training(self, oracle=None):
        if oracle is None:
            oracle = OracleLstm(num_vocabulary=self.vocab_size, batch_size=self.batch_size, emb_dim=self.emb_dim,
                                hidden_dim=self.hidden_dim, sequence_length=self.sequence_length,
                                start_token=self.start_token)
        self.set_oracle(oracle)

        generator = Generator(num_vocabulary=self.vocab_size, batch_size=self.batch_size, emb_dim=self.emb_dim,
                              hidden_dim=self.hidden_dim, sequence_length=self.sequence_length,
                              start_token=self.start_token)
        self.set_generator(generator)

        discriminator = Discriminator(sequence_length=self.sequence_length, num_classes=2, vocab_size=self.vocab_size,
                                      emd_dim=self.emb_dim, filter_sizes=self.filter_size, num_filters=self.num_filters,
                                      l2_reg_lambda=self.l2_reg_lambda)
        self.set_discriminator(discriminator)

        gen_dataloader = DataLoader(batch_size=self.batch_size, seq_length=self.sequence_length)
        oracle_dataloader = DataLoader(batch_size=self.batch_size, seq_length=self.sequence_length)
        dis_dataloader = DisDataloader(batch_size=self.batch_size, seq_length=self.sequence_length)

        self.set_data_loader(gen_loader=gen_dataloader, dis_loader=dis_dataloader, oracle_loader=oracle_dataloader)

    def init_real_training(self, data_loc=None):
        from utils.text_process import text_precess, text_to_code
        from utils.text_process import get_tokenized, get_word_list, get_dict
        if data_loc is None:
            print(f"Error: Training data not specified.")
            sys.exit(1)
        transform_file(data_loc)
        self.sequence_length, self.vocab_size = text_precess(data_loc)
        generator = Generator(num_vocabulary=self.vocab_size, batch_size=self.batch_size, emb_dim=self.emb_dim,
                              hidden_dim=self.hidden_dim, sequence_length=self.sequence_length,
                              start_token=self.start_token)
        self.set_generator(generator)

        discriminator = Discriminator(sequence_length=self.sequence_length, num_classes=2, vocab_size=self.vocab_size,
                                      emd_dim=self.emb_dim, filter_sizes=self.filter_size, num_filters=self.num_filters,
                                      l2_reg_lambda=self.l2_reg_lambda)
        self.set_discriminator(discriminator)

        gen_dataloader = DataLoader(batch_size=self.batch_size, seq_length=self.sequence_length)
        oracle_dataloader = None
        dis_dataloader = DisDataloader(batch_size=self.batch_size, seq_length=self.sequence_length)

        self.set_data_loader(gen_loader=gen_dataloader, dis_loader=dis_dataloader, oracle_loader=oracle_dataloader)
        tokens = get_tokenized(data_loc)
        word_set = get_word_list(tokens)
        [word_index_dict, index_word_dict] = get_dict(word_set)

        with open('save/iw_dict.json', 'w') as file:
            json.dump(index_word_dict, file)

        with open(self.oracle_file, 'w') as outfile:
            outfile.write(text_to_code(tokens, word_index_dict, self.sequence_length))
        return word_index_dict, index_word_dict

    def init_real_metric(self):
        from utils.metrics.DocEmbSim import DocEmbSim
        docsim = DocEmbSim(oracle_file=self.oracle_file, generator_file=self.generator_file,
                           num_vocabulary=self.vocab_size)
        self.add_metric(docsim)

        inll = Nll(data_loader=self.gen_data_loader, rnn=self.generator, sess=self.sess)
        inll.set_name('nll-test')
        self.add_metric(inll)

    def train_real(self, data_loc, model_loc, output_path):
        from utils.text_process import code_to_text
        from utils.text_process import get_tokenized
        wi_dict, iw_dict = self.init_real_training(data_loc)
        self.init_real_metric()

        def get_real_test_file(dict=iw_dict):
            with open(self.generator_file, 'r') as f:
                codes = get_tokenized(self.generator_file)
            output_dir = os.path.dirname(self.output_path)

            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print(output_dir, 'created.')
            with open(self.output_path, 'w') as outfile:
                outfile.write(code_to_text(codes=codes, dictionary=dict))

        self.sess.run(tf.global_variables_initializer())

        self.pre_epoch_num = 80
        self.adversarial_epoch_num = 100

        self.log = open('experiment-log-seqgan-real.csv', 'w')
        generate_samples(self.sess, self.generator, self.batch_size, self.generate_num, self.generator_file)
        self.gen_data_loader.create_batches(self.oracle_file)

        print('start pre-train generator:')
        for epoch in range(self.pre_epoch_num):
            start = time()
            loss = pre_train_epoch(self.sess, self.generator, self.gen_data_loader)
            end = time()
            print('epoch:' + str(self.epoch) + '\t time:' + str(end - start))
            self.add_epoch()
            if epoch % 5 == 0:
                generate_samples(self.sess, self.generator, self.batch_size, self.generate_num, self.generator_file)
                get_real_test_file()
                self.evaluate()

        print('start pre-train discriminator:')
        self.reset_epoch()
        for epoch in range(self.pre_epoch_num):
            print('epoch:' + str(epoch))
            self.train_discriminator()

        self.reset_epoch()
        print('adversarial training:')
        self.reward = Reward(self.generator, .8)
        for epoch in range(self.adversarial_epoch_num):
            start = time()
            for index in range(1):
                samples = self.generator.generate(self.sess)
                rewards = self.reward.get_reward(self.sess, samples, 16, self.discriminator, model_loc)
                feed = {
                    self.generator.x: samples,
                    self.generator.rewards: rewards,
                }
                loss, _ = self.sess.run([self.generator.g_loss, self.generator.g_updates], feed_dict=feed)
                print(loss)
            end = time()
            self.add_epoch()
            print('epoch:' + str(self.epoch) + '\t time:' + str(end - start))
            if epoch % 5 == 0 or epoch == self.adversarial_epoch_num - 1:
                generate_samples(self.sess, self.generator, self.batch_size, self.generate_num, self.generator_file)
                get_real_test_file()
                self.evaluate()

            self.reward.update_params()
            for _ in range(15):
                self.train_discriminator()
