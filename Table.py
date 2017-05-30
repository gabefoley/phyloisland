from flask_table import Table, Col, LinkCol
from flask import url_for

class SortableTable(Table):
    id = Col('ID')
    name = Col('Name')
    species = Col('Species')

    # link = LinkCol(
    #     'Link', 'flask_link', url_kwargs=dict(id='id'), allow_sort=False)
    allow_sort = True

    def sort_url(self, col_key, reverse=False):
        if reverse:
            direction = 'desc'
        else:
            direction = 'asc'
        return url_for('index', sort=col_key, direction=direction)

class Item(object):
      """ a little fake database """
      def __init__(self, id, name, species):
          self.id = id
          self.name = name
          self.species = species


      @classmethod
      def get_elements(cls):
          return itemlist

      def add_elements(cls, item):
          itemlist.append(item)


      @classmethod
      def get_sorted_by(cls, sort, reverse=False):
          return sorted(
              cls.get_elements(),
              key=lambda x: getattr(x, sort),
              reverse=reverse)

      @classmethod
      def get_element_by_id(cls, id):
          return [i for i in cls.get_elements() if i.id == id][0]

itemlist = []
# from flask_table import Table, Col, LinkCol
# from flask import url_for
#
# class SortableTable(Table):
#     id = Col('ID')
#     name = Col('Name')
#     description = Col('Description')
#     link = LinkCol(
#         'Link', 'flask_link', url_kwargs=dict(id='id'), allow_sort=False)
#     allow_sort = True
#
#     def sort_url(self, col_key, reverse=False):
#         if reverse:
#             direction = 'desc'
#         else:
#             direction = 'asc'
#         return url_for('index', sort=col_key, direction=direction)
#
# class Item(object):
#     """ a little fake database """
#     def __init__(self, id, name, description):
#         self.id = id
#         self.name = name
#         self.description = description
#
#     @classmethod
#     def get_elements(cls):
#         return [
#             Item(1, 'Z', 'zzzzz'),
#             Item(2, 'K', 'aaaaa'),
#             Item(3, 'B', 'bbbbb')]
#
#     @classmethod
#     def get_sorted_by(cls, sort, reverse=False):
#         return sorted(
#             cls.get_elements(),
#             key=lambda x: getattr(x, sort),
#             reverse=reverse)
#
#     @classmethod
#     def get_element_by_id(cls, id):
#         return [i for i in cls.get_elements() if i.id == id][0]
