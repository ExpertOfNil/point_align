#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifndef HASH_ALGORITHM
#define HASH_ALGORITHM DJB2
#endif

static inline uint32_t hash_djb2(const char *key, uint32_t capacity) {
    uint32_t hash = 5381;
    int32_t c;
    while ((c = *key++)) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash % capacity;
}
#if HASH_ALGORITHM == DJB2
#define hash(key, capacity) hash_djb2(key, capacity)
#endif

#define DEFAULT_HASHMAP_CAPACITY 16
#define DEFINE_HASHMAP(type)                                                 \
    typedef struct type##Entry type##Entry;                                  \
    struct type##Entry {                                                     \
        char *key;                                                           \
        type value;                                                          \
        type##Entry *next;                                                   \
    };                                                                       \
                                                                             \
    typedef struct type##Map {                                               \
        type##Entry **entries;                                               \
        uint32_t count;                                                      \
        uint32_t capacity;                                                   \
    } type##Map;                                                             \
                                                                             \
    static inline void type##Map_init(type##Map *map) {                      \
        (map)->entries = NULL;                                               \
        (map)->count = 0;                                                    \
        (map)->capacity = DEFAULT_HASHMAP_CAPACITY;                          \
        (map)->entries = calloc((map)->capacity, sizeof(type##Entry *));     \
    }                                                                        \
                                                                             \
    static inline void type##Map_resize(type##Map *map) {                    \
        type##Entry **prev_entries = (map)->entries;                         \
        uint32_t prev_cap = (map)->capacity;                                 \
        (map)->capacity *= 2;                                                \
        (map)->entries = calloc((map)->capacity, sizeof(type##Entry *));     \
        (map)->count = 0;                                                    \
        for (uint32_t i = 0; i < prev_cap; ++i) {                            \
            type##Entry *entry = prev_entries[i];                            \
            while (entry) {                                                  \
                type##Entry *next = entry->next;                             \
                uint32_t new_index = hash_djb2(entry->key, (map)->capacity); \
                entry->next = (map)->entries[new_index];                     \
                (map)->entries[new_index] = entry;                           \
                (map)->count++;                                              \
                entry = next;                                                \
            }                                                                \
        }                                                                    \
        free(prev_entries);                                                  \
    }                                                                        \
                                                                             \
    static inline void type##Map_put(                                        \
        type##Map *map, const char *key, type value                          \
    ) {                                                                      \
        uint32_t index = hash_djb2(key, (map)->capacity);                    \
        type##Entry *entry = (map)->entries[index];                          \
        while (entry) {                                                      \
            if (strcmp(entry->key, key) == 0) {                              \
                entry->value = value;                                        \
                return;                                                      \
            }                                                                \
            entry = entry->next;                                             \
        }                                                                    \
                                                                             \
        if (map->count >= map->capacity) {                                   \
            type##Map_resize(map);                                           \
            index = hash_djb2(key, (map)->capacity);                         \
        }                                                                    \
        type##Entry *new_entry = malloc(sizeof(type##Entry));                \
        new_entry->key = malloc(strlen(key) + 1);                            \
        strcpy(new_entry->key, key);                                         \
        new_entry->value = value;                                            \
        new_entry->next = (map)->entries[index];                             \
        (map)->entries[index] = new_entry;                                   \
        (map)->count++;                                                      \
    }                                                                        \
                                                                             \
    static inline type *type##Map_get(type##Map *map, const char *key) {     \
        uint32_t index = hash_djb2(key, (map)->capacity);                    \
        type##Entry *entry = (map)->entries[index];                          \
        while (entry) {                                                      \
            if (strcmp(entry->key, key) == 0) {                              \
                return &(entry->value);                                      \
            }                                                                \
            entry = entry->next;                                             \
        }                                                                    \
        return NULL;                                                         \
    }                                                                        \
                                                                             \
    static inline void type##Map_clear(type##Map *map) {                     \
        for (uint32_t i = 0; i < (map)->capacity; ++i) {                     \
            type##Entry *entry = (map)->entries[i];                          \
            while (entry) {                                                  \
                type##Entry *next = entry->next;                             \
                free(entry->key);                                            \
                free(entry);                                                 \
                entry = next;                                                \
            }                                                                \
            (map)->entries[i] = NULL;                                        \
        }                                                                    \
        (map)->count = 0;                                                    \
    }                                                                        \
                                                                             \
    static inline void type##Map_free(type##Map *map) {                      \
        for (uint32_t i = 0; i < (map)->capacity; ++i) {                     \
            type##Entry *entry = (map)->entries[i];                          \
            while (entry) {                                                  \
                type##Entry *next = entry->next;                             \
                free(entry->key);                                            \
                free(entry);                                                 \
                entry = next;                                                \
            }                                                                \
        }                                                                    \
        free((map)->entries);                                                \
        (map)->entries = NULL;                                               \
        (map)->count = 0;                                                    \
        (map)->capacity = 0;                                                 \
    }

#endif /* HASHMAP_H */
