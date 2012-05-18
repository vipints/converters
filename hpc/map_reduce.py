import threading
import Queue
import operator
import re

class MapReduce:
    """ MapReduce - to use, subclass by defining these functions,
                    then call self.map_reduce():
        parse_fn(self, k, v) => [(k, v), ...]
        map_fn(self, k, v) => [(k, v1), (k, v2), ...]
        reduce_fn(self, k, [v1, v2, ...]) => [(k, v)]
        output_fn(self, [(k, v), ...])
    """
    def __init__(self):
        self.data = None
        self.num_worker_threads = 5

    class SynchronizedDict(): # we need this for merging
        def __init__(self):
            self.lock = threading.Lock()
            self.d = {}
        def isin(self, k):
            with self.lock:
                if k in self.d:
                    return True
                else:
                    return False
        def get(self, k):
            with self.lock:
                return self.d[k]
        def set(self, k, v): # we don't need del
            with self.lock:
                self.d[k] = v
        def set_append(self, k, v): # for thread-safe list append
            with self.lock:
                self.d[k].append(v)
        def items(self):
            with self.lock:
                return self.d.items() # 

    def create_queue(self, input_list): # helper fn for queues
        output_queue = Queue.Queue()
        for value in input_list:
            output_queue.put(value)
        return output_queue

    def create_list(self, input_queue): # helper fn for queues
        output_list = []
        while not input_queue.empty():
            item = input_queue.get()
            output_list.append(item)
            input_queue.task_done()
        return output_list

    def merge_fn(self, k, v, merge_dict): # helper fn for merge
        if merge_dict.isin(k):
            merge_dict.set_append(k, v)
        else:
            merge_dict.set(k, [v])

    def process_queue(self, input_queue, fn_selector): # helper fn
        output_queue = Queue.Queue()
        if fn_selector == 'merge':
            merge_dict = self.SynchronizedDict()
        def worker():
            while not input_queue.empty():
                (k, v) = input_queue.get()
                if fn_selector in ['map', 'reduce']:
                    if fn_selector == 'map':
                        result_list = self.map_fn(k, v)
                    elif fn_selector == 'reduce':
                        result_list = self.reduce_fn(k, v)
                    for result_tuple in result_list: # flatten
                        output_queue.put(result_tuple)
                elif fn_selector == 'merge': # merge v to same k
                    self.merge_fn(k, v, merge_dict)
                else:
                    raise Exception, "Bad fn_selector="+fn_selector
                input_queue.task_done()
        
        for i in range(self.num_worker_threads): # start threads
            worker_thread = threading.Thread(target=worker)
            worker_thread.daemon = True
            worker_thread.start()
        input_queue.join() # wait for worker threads to finish
        if fn_selector == 'merge':
            output_list = sorted(merge_dict.items(), key=operator.itemgetter(0))
            output_queue = self.create_queue(output_list)
        return output_queue

    def map_reduce(self): # the actual map-reduce algoritm
        data_dict = self.parse_fn(self.data)
        data_queue = self.create_queue(data_dict) # enqueue the data so we can multi-process
        map_queue = self.process_queue(data_queue, 'map') # [(k,v),...] => [(k,v1),(k,v2),...]
        merge_queue = self.process_queue(map_queue, 'merge') # [(k,v1),(k,v2),...] => [(k,[v1,v2,...]),...]
        reduce_queue = self.process_queue(merge_queue, 'reduce') # [(k,[v1,v2,...]),...] => [(k,v),...]
        output_list = self.create_list(reduce_queue) # deque into list for output handling
        self.output_fn(output_list)

"""class WordCount(MapReduce):

    def __init__(self):
        MapReduce.__init__(self)
        self.min_count = 1

    def parse_fn(self, data): # break string into [(k, v), ...] tuples for each line
        data_list = map(lambda line: (None, line), data.splitlines())
        return data_list

    def map_fn(self, key, str): # return (word, 1) tuples for each word, ignore key
        word_list = []
        for word in re.split(r'\W+', str.lower()):
            bare_word = re.sub(r"[^A-Za-z0-9]*", r"", word);
            if len(bare_word) > 0:
                word_list.append((bare_word, 1))
        return word_list

    def reduce_fn(self, word, count_list): # just sum the counts
        return [(word, sum(count_list))]

    def output_fn(self, output_list): # just print the resulting list
        print "Word".ljust(15), "Count".rjust(5)
        print "______________".ljust(15), "_____".rjust(5)
        sorted_list = sorted(output_list, key=operator.itemgetter(1), reverse=True)
        for (word, count) in sorted_list:
            if count > self.min_count:
                print word.ljust(15), repr(count).rjust(5)
        print

    def test_with_monty(self):
        #self.data = The Meaning of Life is:
        #    try and be nice to people,
        #    avoid eating fat,
        #    read a good book every now and then,
        #    get some walking in,
        #    and try and live together in peace and harmony
        #    with people of all creeds and nations.
        self.map_reduce()

class DistributedGrep(MapReduce):

    def __init__(self):
        MapReduce.__init__(self)
        self.matcher = None

    def parse_fn(self, data): # one list item per line with line number
        data_list = []
        line_num = 1
        for line in data.splitlines():
            data_list.append((line_num, line))
            line_num = line_num + 1
        return data_list

    def map_fn(self, line_num, line): # return line if matches, include line num
        matcher = self.matcher
        matched_line = []
        if matcher.match(line):
            matched_line = [(line_num, line)]
        return matched_line

    def reduce_fn(self, line_num, line_list): # identity reducer
        return [(line_num, line_list[0])] # we only ever have one line in the list

    def output_fn(self, output_list): # just print the resulting list
        print "LineNum".rjust(8), "Line".ljust(70)
        print "_______".rjust(8), "____"
        for (line_num, line) in sorted(output_list, key=operator.itemgetter(0)):
            print repr(line_num).rjust(8), line.ljust(70)
        print

def main():
    wc = WordCount()
    wc.test_with_monty()

    #dg = DistributedGrep()

if __name__ == "__main__":
    main()"""

