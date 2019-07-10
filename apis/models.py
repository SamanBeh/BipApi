from django.db import models

class Peptide(models.Model):
    sequence = models.TextField()
    def __str__(self):
        return sequence


class Car(models.Model):
    name = models.CharField(max_length=100)
    top_speed = models.IntegerField()
