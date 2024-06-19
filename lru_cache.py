from collections import OrderedDict


class LRUCache:
    """
    Gotten from https://www.romaglushko.com/blog/design-lru-cache/
    """
    capacity: int
    cache_map: OrderedDict
 
    def __init__(self, capacity: int):
        self.capacity = capacity
        
        self.hits = 0
        self.misses = 0
        self.cache_map = OrderedDict()
 
    def get(self, key: int) -> int:
        try:
            value = self.cache_map[key]
            self.cache_map.move_to_end(key)
            self.hits += 1
            return value
        except KeyError:
            self.misses += 1
            return -1

 
    def put(self, key: int, value: int) -> None:
        if key in self.cache_map:
            self.cache_map[key] = value
            self.cache_map.move_to_end(key)
            return 1

        if len(self.cache_map) >= self.capacity:
            lru_key = next(iter(self.cache_map))
            del self.cache_map[lru_key]

        self.cache_map[key] = value
        return -1

    def size(self):
        return len(self.cache_map)
    
    def cache_info(self):
        return f"Cache info: hits={self.hits}, current_size={self.size()}, capacity={self.capacity}, misses={self.misses}"
    
    def clear(self):
        self.hits = self.misses = 0
        self.cache_map.clear()
