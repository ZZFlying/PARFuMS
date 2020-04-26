def Singleton(cls):
    _instance = {}

    def _singleton(*args, **kargs):
        if cls not in _instance:
            _instance[cls] = cls(*args, **kargs)
        return _instance[cls]

    return _singleton


def read_config(config_file):
    contig = dict()
    # 读取配置文件，忽略#开头的行
    with open(config_file) as config_in:
        for line in config_in:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            array = line.split(':')
            key = array[0].strip()
            value = array[1].strip()
            if value.isdigit():
                value = int(value)
            elif value.upper() == 'TRUE':
                value = True
            elif value.upper() == 'FALSE':
                value = False
            contig[key] = value
    return contig

# 装饰器实现单例模式
@Singleton
class Config(object):
    def __init__(self, config_file=None):
        self._config = read_config(config_file)

    def __getitem__(self, item):
        return self._config[item]

    def __setitem__(self, key, value):
        self._config[key] = value

    def update(self, configs):
        self._config.update(configs)